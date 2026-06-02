# ReactionProfileHunter S4 特征提取完整报告

**版本**: V6.2  
**日期**: 2026-03-16  
**项目**: ReactionProfileHunter (~50k行, 198个Python文件)

---

## 一、S4特征提取架构总览

### 1.1 管道位置

S4 (Step 4) 位于ReactionProfileHunter四步流水线的最后阶段:

```
S1: 构象搜索/锚定 → S2: 逆扫描/前向扫描 → S3: TS优化 → S4: 特征提取
```

**S4的核心职责**: 从S1/S2/S3的量子化学输出中提取反应机理的几何、电子、热化学特征，输出标准化CSV供机器学习使用。

### 1.2 插件化架构

S4采用**插件式特征提取器**架构，核心组件包括:

| 组件 | 文件位置 | 职责 |
|------|---------|------|
| `FeatureMiner` | `feature_miner.py:31` | 主编排器: 发现插件、构建上下文、执行提取、输出3文件 |
| `BaseExtractor` | `extractors/base.py:40` | 插件基类: 定义`extract()`接口和注册机制 |
| `FeatureContext` | `context.py` | 运行时上下文: 路径、QC结果、元数据 |
| `PluginTrace` | `context.py` | 执行轨迹: 状态、警告、错误、耗时 |
| `FeatureSchema` | `schema.py` | Schema管理: FIXED_COLUMNS、MLR_COLUMNS、CSV写入 |

**插件发现机制**:
```python
# extractors/__init__.py 导入所有模块, 触发 register_extractor()
# FeatureMiner 通过 list_extractors() 发现所有已注册插件
EXTRACTORS: Dict[str, BaseExtractor] = {}

def register_extractor(extractor: BaseExtractor) -> None:
    name = extractor.get_plugin_name()
    EXTRACTORS[name] = extractor
```

### 1.3 输出合同 (3文件制)

无论提取成功与否,S4**始终输出**以下3个文件:

```
S4_Data/
├── features_raw.csv      # 完整特征表 (FIXED_COLUMNS + 动态列)
├── features_mlr.csv     # MLR就绪特征表 (默认22列去重版)
└── feature_meta.json     # 元数据: schema版本、配置快照、插件轨迹
```

---

## 二、已注册提取器清单 (V6.2)

| 提取器 | 特征前缀 | 域 | 优先级 |
|--------|---------|-----|--------|
| `thermo` | `thermo.*` | 热力学能 | P0 (核心) |
| `geometry` | `geom.*` | 几何结构 | P0 (核心) |
| `qc_checks` | `qc.*` | QC验证 | P0 (核心) |
| `ts_quality` | `ts.*` | TS质量 | P0 (核心) |
| `step1_activation` | `s1_*` | S1活化特征 | P0 (V6.2新) |
| `step2_cyclization` | `s2_*` | S2环化特征 | P0 (V6.2新) |
| `multiwfn_features` | `mw_*` | 局部反应性 | P2 (Tier-1) |
| `nbo_e2` | `nbo.*` | NBO相互作用 | 可选 |
| `nics` | `arom.*` | 磁性芳香性 | 可选 |
| `interaction_analysis` | `eda.*` | 片段EDA | Phase B |
| `fmo_cdft_dipolar` | `fmo_cdft_dipolar.*` | FMO/CDFT | 默认禁用 |

---

## 三、特征分类与物理/化学意义

### 3.1 热力学特征 (`thermo.*`)

**数据源**: SPMatrixReport / Shermo

| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `thermo.dG_activation` | kcal/mol | **活化自由能** — 反应速率的关键驱动力, k ∝ exp(-dG‡/RT) |
| `thermo.dG_reaction` | kcal/mol | **反应自由能** — 反应放热性/热力学驱动力 |
| `thermo.dE_activation` | kcal/mol | 电子活化能 ΔE(TS) - ΔE(R) |
| `thermo.dE_reaction` | kcal/mol | 电子反应能 ΔE(P) - ΔE(R) |
| `thermo.energy_source` | - | 能量来源标签 ("gibbs" 或 "electronic") |

**化学意义**: 活化能决定反应速率, 反应能决定热力学可行性。Gibbs能优先(包含熵贡献), 无Gibbs时降级到电子能。

### 3.2 几何特征 (`geom.*`)

**数据源**: TS XYZ坐标 + forming_bonds

| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `geom.r1`, `geom.r2` | Å | **形成键长度** — 过渡态中正在形成/断裂的键 |
| `geom.asynch` | Å | **不对称性** = \|r1 - r2\| — 键形成的同步程度 |
| `geom.asynch_index` | - | 归一化不对称性 = asynch / (r1 + r2) |
| `geom.r_avg` | Å | 平均形成键长 — 键形成的整体进程 |
| `geom.dr` | 有符号 | **有符号不对称性** = r1 - r2 — 反映键形成顺序 |
| `geom.rg_ts` | Å | **转动半径** — TS结构紧凑程度 |
| `geom.close_contacts` | - | 近距接触数 (< 2.2 Å) — 立体压力指标 |
| `geom.close_contacts_density` | - | 接触密度 = close_contacts / natoms — 单位原子立体压力 |

**化学意义**: 过渡态几何反映反应机理 (协同 vs 分步)、立体电子效应、环张力等。asynch是区分 concerted/sequential 机理的关键指标。

### 3.3 Step 1 激活特征 (`s1_*`)

**数据源**: S1输出 (Shermo/HOAc热化学)

| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `s1_dG_act` | kcal/mol | Step 1活化吉布斯能 — 底物→活化中间体 |
| `s1_Keq_act` | - | 平衡常数 = exp(-dG/RT) |
| `s1_Nconf_eff` | - | **有效构象数** (Boltzmann加权) — 构象多态性 |
| `s1_Sconf` | cal/mol·K | **构象熵** — 构象灵活性 |
| `s1_E_avg_weighted` | kcal/mol | Boltzmann加权平均能量 |
| `s1_E_std` | kcal/mol | 能量标准差 — 构象能隙 |
| `s1_tau_CH_C_O` | deg | 离去基团二面角 (C-H···C-O) — α-H门控效应 |

**化学意义**: 底物活化障碍和构象多态性影响反应活性和选择性。构象熵在低温反应中尤为重要。

### 3.4 Step 2 环化特征 (`s2_*`)

**数据源**: S3 TS相关文件 (fchk/log/ORCA输出)

#### 3.4.1 动力学/热化学
| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `s2_dGddagger` | kcal/mol | 环化TS活化吉布斯能 |
| `s2_dHddagger` | kcal/mol | 活化焓 |
| `s2_dSddagger` | cal/mol·K | 活化熵 |
| `s2_TdSddagger` | kcal/mol | TΔS‡ 温度修正项 |
| `s2_dGrxn` | kcal/mol | 环化反应吉布斯能 |

#### 3.4.2 CDFT (概念DFT) 特征
| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `s2_eps_homo` | eV | **HOMO能量** — 给电子倾向 (亲核性) |
| `s2_eps_lumo` | eV | **LUMO能量** — 受电子倾向 (亲电性) |
| `s2_mu` | eV | **化学势** = (HOMO+LUMO)/2 — 电子流向 |
| `s2_eta` | eV | **硬度** = (LUMO-HOMO)/2 — 电子供给抗性 |
| `s2_omega` | eV | **电吸电性** = μ²/2η — 受电子能力 |

**化学意义**: CDFT提供反应物的整体电子结构描述, 预测极性反应的方向和驱动力。omega (electrophilicity) 是筛选亲电试剂的关键指标。

#### 3.4.3 GEDT (全局电子密度转移)
| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `s2_gedt_value` | e | 从Nucleophile→Electrophile转移的电子数 |

**化学意义**: GEDT量化 cycloaddition 中的电荷转移程度, 区分 normal electron demand (NED) vs inverse electron demand (IED) 反应。

#### 3.4.4 TS有效性
| 特征名 | 化学意义 |
|--------|----------|
| `s2_ts_validity_flag` | 0/1 — TS质量检查 (恰好1个虚频) |

### 3.5 TS质量特征 (`ts.*`)

**数据源**: TS Gaussian log / ORCA输出

| 特征名 | 单位 | 化学意义 |
|--------|------|----------|
| `ts.n_imag` | - | **虚频个数** (理想值=1) |
| `ts.imag1_cm1_abs` | cm⁻¹ | **最负虚频绝对值** — 反应坐标的清晰度 |
| `ts.dipole_debye` | Debye | **偶极矩** — TS极性 |

**化学意义**: 虚频是TS的确认标准; 虚频大小反映反应坐标的曲率(陡峭 vs 平坦)。

### 3.6 Multiwfn特征 (`mw_*`)

**数据源**: Multiwfn分析 (fchk文件)

| 特征名 | 化学意义 |
|--------|----------|
| `mw_fukui_f+_forming1` | **Fukui阳电子体指数** — 亲核性 |
| `mw_fukui_f-_forming1` | **Fukui阴电子体指数** — 亲电性 |
| `mw_fukui_f0` | 中性指数 = (f+ + f-)/2 |
| `mw_dual_descriptor` | **双描述符** = f+ - f- — 化学活性 |
| `mw_rho_bcp` | QTAIM电子密度 — 键强度 |

**化学意义**: Fukui函数提供**局部**反应性描述, 识别反应位点。与CDFT的**整体**描述互补。

### 3.7 QC验证特征 (`qc.*`)

| 特征名 | 化学意义 |
|--------|----------|
| `qc.has_gibbs` | 是否有Gibbs能 |
| `qc.used_fallback_electronic` | 是否降级到电子能 |
| `qc.sp_report_validated` | SPMatrixReport验证通过 |
| `qc.forming_bonds_valid` | forming_bonds有效 |
| `qc.warnings_count` | 全局警告计数 |
| `qc.sample_weight` | **样本权重** (1.0=严格QC, 0.0=降级) |

---

## 四、MLR就绪特征集 (降维策略)

### 4.1 默认22列 (V6.2去重版)

```python
MLR_COLUMNS_V3_DEDUP = [
    "sample_id",
    # 热力学
    "thermo.dE_activation",
    "thermo.dE_reaction",
    # 几何
    "geom.r_avg",
    "geom.dr",
    "geom.close_contacts_density",
    # TS质量
    "ts.n_imag",
    "ts.imag1_cm1_abs",
    # Step1
    "s1_dG_act",
    "s1_Keq_act",
    "s1_Nconf_eff",
    "s1_Sconf",
    "s1_E_avg_weighted",
    "s1_E_std",
    "s1_tau_CH_C_O",
    # Step2电子结构
    "s2_eps_homo",
    "s2_eps_lumo",
    "s2_mu",
    "s2_eta",
    "s2_omega",
    "s2_gedt_value",
    # 样本质量
    "qc.sample_weight",
]
```

### 4.2 特征维度分析

| 类别 | 原始特征数 | MLR保留数 | 去除原因 |
|------|-----------|----------|----------|
| thermo | 10 | 2 | 冗余 (dG/dE, Gibbs/Electronic) |
| geom | 13 | 3 | 冗余 (r1/r2 → r_avg, dr) |
| ts | 3 | 2 | 保留关键 |
| s1 | 7 | 7 | 全保留 |
| s2 | 20+ | 7 | 仅保留电子结构, 去重几何 |
| qc | 6 | 1 | 仅保留权重 |
| **总计** | **~60** | **22** | **去重63%** |

---

## 五、机器学习重要性与维度灾难

### 5.1 维度灾难问题

S4当前提供约**60个原始特征**, 经过去重后MLR集有**22个特征**。这引发以下问题:

1. **过拟合风险**: 22维特征空间 vs 典型数据集大小 (100-1000个样本)
2. **特征共线性**: 多个特征高度相关 (如 `s2_eps_homo` + `s2_eps_lumo` → `s2_eta`, `s2_mu`)
3. **噪声放大**: 低质量特征 (NaN, 高方差) 污染模型

### 5.2 当前缓解策略

1. **特征去重** (V6.2): 移除与FIXED_COLUMNS重复的s2几何特征
2. **样本权重** (`qc.sample_weight`): 标记低质量样本, 训练时可降权
3. **降级处理**: 缺失特征填充NaN, 保留信息但允许模型处理

### 5.3 建议的进一步降维

| 方法 | 描述 | 预期收益 |
|------|------|----------|
| **相关性过滤** | 移除 \|r\| > 0.9 的特征对 | 去除冗余 |
| **PCA/SVD** | 主成分分析, 保留95%方差 | 正交化 |
| **L1正则化** | Lasso回归进行特征选择 | 稀疏解 |
| **树模型重要性** | RF/XGBoost特征重要性排序 | 解释性筛选 |
| **领域知识** | 化学直觉筛选关键描述符 | 最可靠 |

---

## 六、未来机器学习算法建议

### 6.1 回归任务 (预测 dG‡)

| 算法 | 优势 | 劣势 | 适用场景 |
|------|------|------|----------|
| **XGBoost/LightGBM** | 梯度提升, 处理缺失值, 特征重要性 | 需要调参 | 中小数据集, 快速baseline |
| **随机森林** | 鲁棒, 抗过拟合, 特征重要性 | 低方差解释 | 特征筛选初期 |
| **图神经网络 (GNN)** | 分子图直接学习, 端到端 | 需要大数据集 | 新反应类型泛化 |
| **Transformer (MolBERT)** | 预训练, SMILES/3D表示 | 计算成本高 | 大规模数据 |
| **贝叶斯回归** | 不确定性量化, 小数据友好 | 推理慢 | 关键反应预测 |

### 6.2 分类任务 (反应类型/机理)

| 算法 | 目标 |
|------|------|
| **逻辑回归 + L1** | 机理分类 (协同/分步) |
| **SVM** | 小样本高维分类 |
| **GCN/GAT** | 端到端反应分类 |

### 6.3 推荐技术路线

**Phase 1: 快速Baseline**
```
1. 22维MLR特征 → XGBoost回归 → R²评估
2. 相关性过滤 → 特征数降至15以内
3. 样本权重应用 → 评估低质量样本影响
```

**Phase 2: 增强模型**
```
1. 加入mw_* (Fukui) 特征 → 扩展到30维
2. 尝试 GNN (GCN/GAT) → 分子图表示
3. 集成学习 → XGB + RF + NN
```

**Phase 3: 深度学习**
```
1. 3D分子构型 → SchNet/Equivariant Transformer
2. 反应模板 → RxnFP fingerprint
3. 大规模预训练 → ChemBERTa/MolBERT
```

### 6.4 关键注意事项

1. **数据增强**: 构象平均、能量扰动 → 扩充训练集
2. **不确定性**: 贝叶斯方法或MC Dropout → 预测置信度
3. **可解释性**: SHAP/LIME → 特征贡献可视化
4. **域外检测**: 预测分布外样本 → 拒绝服务

---

## 七、附录: 关键文件索引

| 类别 | 文件路径 |
|------|----------|
| 主编排 | `rph_core/steps/step4_features/feature_miner.py` |
| 插件基类 | `rph_core/steps/step4_features/extractors/base.py` |
| Schema定义 | `rph_core/steps/step4_features/schema.py` |
| 上下文 | `rph_core/steps/step4_features/context.py` |
| 热力学提取器 | `rph_core/steps/step4_features/extractors/thermo.py` |
| 几何提取器 | `rph_core/steps/step4_features/extractors/geometry.py` |
| Step1提取器 | `rph_core/steps/step4_features/extractors/step1_activation.py` |
| Step2提取器 | `rph_core/steps/step4_features/extractors/step2_cyclization.py` |
| Multiwfn | `rph_core/steps/step4_features/extractors/multiwfn_features.py` |
| 配置 | `config/defaults.yaml` (step4.mlr.columns) |

---

**文档版本**: V1.0  
**完成时间**: 2026-03-16  
**数据来源**: 代码分析 + S4_FEATURES_SUMMARY.md + schema.py
