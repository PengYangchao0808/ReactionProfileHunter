# Step 4: 特征提取总结 (V6.2)

## 概览
ReactionProfileHunter Step 4（S4）是特征提取和打包模块，从S1、S2、S3的结构与QC输出中提取反应机理的几何、电子、热化学特征。

**版本**: V6.2 (extract-only plugin pipeline)  
**输出合同**: 3文件制
- `features_raw.csv` - 完整特征表（所有可用特征）
- `features_mlr.csv` - MLR就绪特征表（固定24列，V6.2）
- `feature_meta.json` - 元数据、schema版本、运行溯源

---

## I. 特征分类与物理意义

### 1. **热力学特征** (thermo.*)

#### 源头: SPMatrixReport / Shermo

| 特征名 | 单位 | 物理意义 | 备注 |
|--------|------|--------|------|
| `thermo.dG_activation` | kcal/mol | **活化自由能** (Gibbs优先，电子能备选) | 反应速率的关键 |
| `thermo.dG_reaction` | kcal/mol | **反应自由能** (Gibbs优先，电子能备选) | 反应放热性 |
| `thermo.dG_activation_gibbs` | kcal/mol | Gibbs活化能 (显式) | 精确热力学 |
| `thermo.dG_reaction_gibbs` | kcal/mol | Gibbs反应能 (显式) | 精确热力学 |
| `thermo.dE_activation` | kcal/mol | **电子活化能** | ΔE(TS) - ΔE(R) |
| `thermo.dE_reaction` | kcal/mol | **电子反应能** | ΔE(P) - ΔE(R) |
| `thermo.energy_source_activation` | - | 能量来源标签 | "gibbs" 或 "electronic" |
| `thermo.energy_source_reaction` | - | 能量来源标签 | "gibbs" 或 "electronic" |
| `thermo.method` | - | 计算方法 (B3LYP/6-31G等) | 跟踪方法 |
| `thermo.solvent` | - | 溶剂标识 | 溶液效应背景 |

**物理意义**: 反应的热力学难度和产热性

---

### 2. **几何特征** (geom.*)

#### 源头: TS XYZ坐标 + forming_bonds

| 特征名 | 单位 | 物理意义 | 备注 |
|--------|------|--------|------|
| `geom.natoms_ts` | - | **TS原子数** | 反应体系规模 |
| `geom.r1` | Å | **第1个形成键长度** | 第1对原子间距 |
| `geom.r2` | Å | **第2个形成键长度** | 第2对原子间距 |
| `geom.asynch` | Å | **不对称性** = \|r1 - r2\| | TS形成键同步程度 |
| `geom.asynch_index` | - | **相对不对称性** = asynch / (r1 + r2) | 归一化不对称性 (0-1) |
| `geom.r_avg` | Å | **平均形成键长** = (r1 + r2) / 2 | 键形成进程 |
| `geom.dr` | Å | **有符号不对称性** = r1 - r2 | 反映键形成顺序 |
| `geom.rg_ts` | Å | **转动半径** | TS结构紧凑程度 |
| `geom.min_nonbonded` | Å | **最小非键距** | 分子堆积程度 |
| `geom.close_contacts` | - | **近距接触数** (< 2.2 Å) | 立体压力 |
| `geom.close_contacts_density` | - | **接触密度** = close_contacts / natoms_ts | 单位原子立体压力 |

**物理意义**: TS的几何结构和立体限制

---

### 3. **Step 1 激活特征** (s1_*) [V6.2新增]

#### 源头: S1输出 (Shermo/HOAc热化学)

| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `s1_dG_act` | kcal/mol | Step 1 **活化吉布斯能** (底物→活化中间体) |
| `s1_Keq_act` | - | Step 1 **平衡常数** = exp(-dG/RT) |
| `s1_Nconf_eff` | - | **有效构象数** = Boltzmann加权 |
| `s1_Sconf` | cal/mol·K | **构象熵** (从构象群体) |
| `s1_E_avg_weighted` | kcal/mol | **Boltzmann加权平均能量** |
| `s1_E_std` | kcal/mol | **加权能量标准差** |
| `s1_tau_CH_C_O` | deg | **离去基团二面角** (C-H···C-O) |

**物理意义**: 底物活化障碍和构象多态性

---

### 4. **Step 2 环化特征** (s2_*) [V6.2新增]

#### 源头: S3 TS相关文件 (fchk/log/ORCA输出)

#### 4a. 动力学特征
| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `s2_dGddagger` | kcal/mol | **TS活化吉布斯能** (环化反应) |
| `s2_dHddagger` | kcal/mol | TS活化焓 |
| `s2_dSddagger` | cal/mol·K | TS活化熵 |
| `s2_TdSddagger` | kcal/mol | 温度修正项 TΔS‡ |
| `s2_dGrxn` | kcal/mol | 环化反应的吉布斯能 |

#### 4b. TS几何特征
| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `s2_d_forming_1` | Å | **第1个形成键距** |
| `s2_d_forming_2` | Å | **第2个形成键距** |
| `s2_asynch` | - | **键形成不对称性指数** |

#### 4c. CDFT (Conceptual DFT) 特征
| 特征名 | 单位 | 物理意义 | 化学含义 |
|--------|------|--------|--------|
| `s2_eps_homo` | eV | **HOMO能量** | 给电子倾向 |
| `s2_eps_lumo` | eV | **LUMO能量** | 受电子倾向 |
| `s2_mu` | eV | **化学势** = (HOMO+LUMO)/2 | 电子流向 |
| `s2_eta` | eV | **硬度** = (LUMO-HOMO)/2 | 电子供给抗性 |
| `s2_omega` | eV | **电吸电性 (Electrophilicity)** = μ²/2η | 受电子能力 |

**CDFT物理意义**: 反应物的极性和电子供给-受体能力

#### 4d. GEDT (全局电子密度转移)
| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `s2_gedt_value` | e | **从Nucleophile→Electrophile转移的电子数** | 荷电转移程度 |

#### 4e. TS有效性
| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `s2_ts_validity_flag` | 0/1 | **TS质量检查** (恰好1个虚频?) |

---

### 5. **TS质量特征** (ts.*)

#### 源头: TS Gaussian log / ORCA输出

| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `ts.n_imag` | - | **虚频个数** (理想值=1) |
| `ts.imag1_cm1_abs` | cm⁻¹ | **最负虚频绝对值** | 反应坐标的清晰度 |
| `ts.dipole_debye` | Debye | **偶极矩大小** | TS的极性 |

**物理意义**: TS的优化质量和化学合理性

---

### 6. **Multiwfn特征** (mw_*) [V6.2可选]

#### 源头: Multiwfn分析 (fchk文件，缓存机制)

| 特征名 | 单位 | 物理意义 |
|--------|------|--------|
| `mw_fukui_f+_forming1` | - | **Fukui阳电子体指数** (原子1) | 亲核性 |
| `mw_fukui_f-_forming1` | - | **Fukui阴电子体指数** (原子1) | 亲电性 |
| `mw_fukui_f0_forming1` | - | **Fukui中性指数** = (f+ + f-)/2 |  |
| `mw_dual_descriptor_forming1` | - | **双描述符** = f+ - f- | 化学活性 |
| `mw_dual_descriptor_forming2` | - | **双描述符** (原子2) | 化学活性 |
| `mw_status` | str | Multiwfn运行状态 | "ok" / "skipped" / "failed" |
| `mw_cache_hit` | bool | 缓存命中 | 性能指标 |

**物理意义**: TS原子的局部反应性 (Fukui函数)

---

### 7. **QC验证特征** (qc.*)

#### 源头: SPMatrixReport + forming_bonds / Feature context

| 特征名 | 范围 | 物理意义 |
|--------|------|--------|
| `qc.has_gibbs` | 0/1 | 是否有Gibbs能 |
| `qc.used_fallback_electronic` | 0/1 | 是否降级到电子能 |
| `qc.sp_report_validated` | 0/1 | SPMatrixReport通过验证 |
| `qc.forming_bonds_valid` | 0/1 | forming_bonds有效 |
| `qc.warnings_count` | int | 全局警告计数 |
| `qc.sample_weight` | 0.0-1.0 | **样本权重** (严格QC检查可降级) |

**物理意义**: 数据质量指标

---

## II. MLR就绪特征 (去重默认集)

默认列集 (MLR_COLUMNS_V3_DEDUP, 22列；也可由 config.step4.mlr.columns 覆盖):

```python
[
    "sample_id",                          # 样本标识
    "thermo.dE_activation",               # 电子活化能
    "thermo.dE_reaction",                 # 电子反应能
    "geom.r_avg",                         # 平均形成键长
    "geom.dr",                            # 有符号不对称性
    "geom.close_contacts_density",        # 接触密度
    "ts.n_imag",                          # 虚频个数
    "ts.imag1_cm1_abs",                   # 最负虚频绝对值
    
    # Step1 特征
    "s1_dG_act",                          # Step1活化能
    "s1_Keq_act",                         # 平衡常数
    "s1_Nconf_eff",                       # 有效构象数
    "s1_Sconf",                           # 构象熵
    "s1_E_avg_weighted",                  # 加权平均能量
    "s1_E_std",                           # 能量标准差
    "s1_tau_CH_C_O",                      # 离去基团二面角
    
    # Step2 特征 (仅保留电子/电荷转移；TS几何与虚频由 geom.* / ts.* 提供)
    "s2_eps_homo",                        # HOMO能量
    "s2_eps_lumo",                        # LUMO能量
    "s2_mu",                              # 化学势
    "s2_eta",                             # 硬度
    "s2_omega",                           # 电吸电性
    "s2_gedt_value",                      # GEDT值

    "qc.sample_weight",                   # 样本权重
]
```

---

## III. 提取器架构 (Plugin Pipeline)

### 已注册提取器 (V6.2)

| 提取器 | 特征前缀 | 类型 | 状态 |
|--------|--------|------|------|
| **thermo** | `thermo.*` | 热力学 | ✓ 核心 |
| **geometry** | `geom.*` | 几何 | ✓ 核心 |
| **qc_checks** | `qc.*` | QC验证 | ✓ 核心 |
| **ts_quality** | `ts.*` | TS质量 | ✓ 核心 |
| **step1_activation** | `s1_*` | 机理感知 | ✓ V6.2新 |
| **step2_cyclization** | `s2_*` | 机理感知 | ✓ V6.2新 |
| **fmo_cdft_dipolar** | `fmo_cdft_dipolar.*` | 电子结构 | ✗ 默认禁用（与 s2_* 重复） |
| **multiwfn_features** | `mw_*` | 局部反应性 | ⚙️ P2 (Tier-1) |
| **interaction** | `eda.*` | 相互作用 | 🔄 Phase B (job_specs) |
| **nbo_e2** | `nbo.*` | NBO | ⚠️ 可选 |
| **nics** | `nics.*` | 磁性 | ⚠️ 可选 |
| **asm_enrichment** | `asm.*` | 活性空间 | ✗ 已取消（默认禁用） |

### 提取器生命周期

```
FeatureContext (路径 + forming_bonds等)
        ↓
  validate_inputs() [检查依赖]
        ↓
     extract() [主逻辑]
        ↓
   PluginTrace [运行轨迹]
  - status: OK/PARTIAL/FAILED/SKIPPED
  - extracted_features: {feature_name: value}
  - warnings, errors, job_specs
        ↓
  FeatureResult (聚合)
        ↓
  outputs:
    - features_raw.csv (所有特征)
    - features_mlr.csv (config 指定列；默认去重 22 列)
    - feature_meta.json (元数据)
```

---

## IV. 输出文件格式

### 4.1 features_raw.csv

```
sample_id, thermo.dG_activation, thermo.dE_activation, ..., s2_omega, mw_fukui_f+_forming1, ...
rx_001,    18.5,                 21.3,                ..., 3.2,      0.128,                ...
rx_002,    16.2,                 19.8,                ..., 2.9,      0.135,                ...
```

- **列顺序**: 固定列 (FIXED_COLUMNS) → 动态列 (sorted)
- **NaN处理**: 缺失特征填充 NaN
- **精度**: float64

### 4.2 features_mlr.csv

```
sample_id, thermo.dE_activation, ..., qc.sample_weight
rx_001,    21.3,                 ..., 1.0
rx_002,    19.8,                 ..., 0.9
```

- **列**: 由 `config.step4.mlr.columns` 指定（defaults.yaml 默认去重 22 列）
- **用途**: 直接输入机器学习管线
- **样本权重**: QC严格性降级标记

### 4.3 feature_meta.json

```json
{
  "schema_version": "6.2",
  "schema_signature": "abc123def456...",
  "feature_status": "ok",
  "method": "extract-only",
  "temperature_K": 298.15,
  "enabled_plugins": ["thermo", "geometry", "qc_checks", "ts_quality", "step1_activation", "step2_cyclization"],
  "plugin_traces": {
    "thermo": {
      "status": "ok",
      "runtime_ms": 12.3,
      "extracted_features": {...},
      "warnings": [],
      "errors": []
    },
    ...
  },
  "config_snapshot": {...},
  "provenance": {
    "extract_mode": "extract-only",
    "plugin_pipeline_version": "6.2",
    "multiwfn_status": "unknown"
  }
}
```

---

## V. 关键物理量总结表

### 反应速率控制

| 量 | 单位 | 意义 | 指导作用 |
|----|------|------|---------|
| dG‡ | kcal/mol | **活化吉布斯能** | 决定反应速率 k ∝ exp(-dG‡/RT) |
| n_imag | - | 虚频数 (理想=1) | TS质量检查 |
| asynch | Å | 键形成不对称性 | 反应机理可视化 |

### 反应放热性

| 量 | 单位 | 意义 |
|----|------|------|
| dG_rxn | kcal/mol | 反应自由能 (>0=不利，<0=有利) |
| dE_rxn | kcal/mol | 反应电子能 |

### 立体因素

| 量 | 单位 | 意义 | 影响 |
|----|------|------|------|
| close_contacts_density | - | 单位体积接触数 | 立体压力 → dG‡↑ |
| rg_ts | Å | TS紧凑度 | 分子大小，热容量 |

### 电子因素

| 量 | 单位 | 意义 |
|----|------|------|
| eta (CDFT) | eV | 轨道间隙 | 反应物稳定性 |
| omega | eV | 电吸电性 | 亲电性 |
| Fukui f+/f- | - | 反应位点 | 可视化反应活性 |

---

## VI. 常见问题与缺失特征处理

### 为什么某个特征是 NaN?

| 情况 | 原因 | 日志记录 |
|------|------|---------|
| `thermo.*` NaN | 缺失SPMatrixReport或g_*/e_*值 | WARNING: 无热力学信息 |
| `geom.r1/r2` NaN | 缺失forming_bonds或TS坐标 | WARNING: 无几何数据 |
| `s1_*` NaN | 缺失S1输出目录 | WARNING: S1 Shermo缺失 |
| `s2_*` NaN | 缺失S3 fchk或log | WARNING: S3 fchk缺失 |
| `mw_*` NaN | Multiwfn禁用或执行失败 | DEBUG: Multiwfn跳过 |
| `s2_ts_validity_flag`=0 | 虚频≠1 或其他QC问题 | WARNING: TS虚频异常 |

### 降级策略 (Graceful Degradation)

- ✓ 缺失Gibbs能 → 降级到电子能 (`energy_source_activation="electronic"`)
- ✓ 缺失S1特征 → `s1_*` = NaN, 继续提取其他
- ✓ Multiwfn失败 → `mw_status="failed"`, 非阻断
- ✓ 样本权重 → `qc.sample_weight` 记录可信度 (1.0=高, 0.0=降级)

---

## VII. 使用示例

### 读取特征

```python
import pandas as pd
import json

# 读取特征表
df_raw = pd.read_csv("S4_Data/features_raw.csv")
df_mlr = pd.read_csv("S4_Data/features_mlr.csv")

# 读取元数据
with open("S4_Data/feature_meta.json") as f:
    meta = json.load(f)

# 检查启用的提取器
print(meta['enabled_plugins'])

# 检查样本质量
print(df_mlr[['sample_id', 's2_ts_validity_flag', 'qc.sample_weight']])
```

### 基本分析

```python
# 活化能分布
import matplotlib.pyplot as plt
plt.hist(df_mlr['thermo.dE_activation'], bins=20)
plt.xlabel('dE‡ (kcal/mol)')
plt.title('活化能分布')

# 相关性: 活化能 vs CDFT电吸电性
plt.scatter(df_mlr['s2_omega'], df_mlr['thermo.dE_activation'])
plt.xlabel('omega (eV)')
plt.ylabel('dE‡ (kcal/mol)')
plt.title('电吸电性 vs 活化能')
```

---

## VIII. 参考资源

- **代码位置**:
  - 主编排: `rph_core/steps/step4_features/feature_miner.py`
  - 提取器实现: `rph_core/steps/step4_features/extractors/*.py`
  - Schema定义: `rph_core/steps/step4_features/schema.py`
  - 上下文: `rph_core/steps/step4_features/context.py`

- **相关文档**:
  - [step4_features/AGENTS.md](rph_core/steps/step4_features/AGENTS.md) - 详细架构
  - [README.md](README.md) - 管线概览
  - [QUICK_START_DUAL_LEVEL.md](QUICK_START_DUAL_LEVEL.md) - 双层策略

---

**文档版本**: V6.2 (2026-02-02)  
**最后更新**: [用户生成的整理文档]  
**维护者**: RPH Team
