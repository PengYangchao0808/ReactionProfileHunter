# Lewis 酸强度拟合与三模式完善方案

**日期**: 2026-07-06  
**状态**: 实施设计稿  
**目标仓库**: ReactionProfileHunter v3.0.1  
**适用范围**: RPH 的 S0-S3 DFT 计算流程。S4 特征提取与 ML 训练仍由 RPH_Postprocess 等外部仓库负责。

## 1. 目标概述

RPH 当前已经具备 Lewis 酸体系的基础实现：可以在几何层面把 Lewis 酸 surrogate 接入反应体系，并在 S2/S3 对配位质量和 TS 虚频质量做门控。下一步需要完成两件事。

第一，将用户侧 Lewis 酸入口收敛为三种清晰模式：

| 用户模式 | 内部建议值 | 是否启用 LA | 默认模型 | 化学含义 |
|---|---|---:|---|---|
| 无 Lewis 酸 | `none` | 否 | 无 | 原始无添加剂通道 |
| 强 Lewis 酸 | `strong_single_site` | 是 | `BF3` | 单配位中心，默认三氟化硼模型 |
| 弱 Lewis 酸 | `weak_multi_site` | 是 | `MgCl2` | 多配位中心，默认氯化镁模型 |

第二，在小分子模块中加入可复用的 Lewis 酸强度拟合流程。该流程以丙酮为标准探针分子，计算丙酮与 Lewis 酸配位前后羰基碳中心电荷的变化：

```text
delta_q_c = q_c(acetone-LA complex) - q_c(free acetone)
```

其中 `delta_q_c` 越正，代表 Lewis 酸使羰基碳更加缺电子，操作性 Lewis 酸强度越高。该指标用于描述真实 Lewis 酸强弱差异，并写入 registry 和 small-molecule cache，供后续外部后处理或 ML 使用。

本方案不把 Lewis 酸强度直接写成反应能垒校正项，也不在 RPH 内新增 S4 特征挖掘逻辑。RPH 负责 DFT 计算、校准结果落盘和元数据传递。

## 2. 当前实现情况

### 2.1 配置层

主配置 `config/defaults.yaml` 已包含 Lewis 酸模块：

```yaml
lewis_acid:
  enabled: false
  mode: "neutral_metal_chloride_chelation_surrogate"
  surrogate: "LiCl"
  surrogate_registry:
    LiCl:
      name: "LiCl"
      smiles: "[Li]Cl"
      charge: 0
      multiplicity: 1
    MgCl2:
      name: "MgCl2"
      smiles: "Cl[Mg]Cl"
      charge: 0
      multiplicity: 1
    AlCl3:
      name: "AlCl3"
      smiles: "Cl[Al](Cl)Cl"
      charge: 0
      multiplicity: 1
    BF3:
      name: "BF3"
      smiles: "B(F)(F)F"
      charge: 0
      multiplicity: 1
```

现有能力：

- 已有 `LiCl`、`MgCl2`、`AlCl3`、`BF3` surrogate registry。
- 已有 `surrogate_overrides`，可以按 surrogate 调整 S2/S3 quality gate。
- `BF3` 和 `MgCl2` 已能表达不同质量门控倾向。

当前不足：

- `mode` 仍是历史内部字符串，不是面向用户的稳定枚举。
- 用户需要直接理解并选择 `surrogate`，没有统一的强/弱/无三档接口。
- 默认 `LiCl` 仍偏向早期中性盐配位探针设计，与当前强酸 `BF3`、弱酸 `MgCl2` 的新目标不完全一致。
- `LiCl`、`AlCl3` 等需要保留为兼容或高级模式，但不应混入常规三档入口。

### 2.2 数据模型

`rph_core/lewis_acid/models.py` 已定义：

- `LewisAcidAdditive`
- `LiCoordinationTrace`
- `LewisAcidQualityFlags`

`LewisAcidAdditive` 当前记录：

- 是否启用
- surrogate 名称与 SMILES
- 电荷与自旋多重度
- additive 原子索引
- additive 元素列表
- organic 原子数
- 中心原子索引，字段名仍为 `metal_atom_index`
- 配置快照

优点：

- 已能把 additive 原子索引从有机分子原子索引中隔离出来。
- 已能避免 S2/S3 把 Lewis 酸原子误当作形成键原子。
- `center_element`、`center_atom_index` 已作为兼容别名存在。

不足：

- `LiCoordinationTrace` 和 `li_o_distances`、`li_cl_distances` 等命名仍带有 LiCl 历史痕迹。
- `metal_cl_max` 对 BF3 并不准确，因为 BF3 是 B-F 中心配体结构，不是 metal-Cl。
- 多配位弱 Lewis 酸只记录最近 O 位点，难以表达多个竞争配位中心。

### 2.3 Orchestrator 接入

`ReactionProfileHunter._prepare_lewis_acid_system()` 当前行为：

1. 读取 `lewis_acid` 配置。
2. 如果未启用，返回原始 product 和 precursor SMILES。
3. 如果启用，将 surrogate SMILES 拼接到几何 SMILES。
4. 保持 canonical product SMILES 不变。
5. 返回 `LewisAcidAdditive`。

S1 结束后，orchestrator 会检查 product XYZ 末尾是否为 additive 原子，并写入：

- `additive_atom_indices`
- `organic_atom_count`
- `metal_atom_index`

这是正确的主线设计，应继续保留。

### 2.4 S2 质量门控

`RetroScanner` 当前实现了 LA 配位轨迹和质量门控：

- 中心原子与卤素或配体距离过大，标记 `la_dissociation`
- 中心原子与底物 O 距离过大，标记 `la_no_coordination`
- 配位位点发生切换且配置不允许，标记 `la_coordination_switch`

这些结果会写入 `LewisAcidQualityFlags`，并可使 S2 结果降级为 `DEGRADED`。

当前不足：

- 日志和字段仍写作 Li-Cl、Li-O，不适合 BF3 和 MgCl2 的统一表达。
- 对 MgCl2 这类弱多配位模型，应允许并记录配位切换，而不是只给布尔值。
- 对 BF3 应重点判断 B-O 配位与 B-F 结构完整性。

### 2.5 S3 质量门控

`TSValidator.validate_la_mode_projection()` 当前判断 TS 虚频是否真正落在形成键反应坐标上：

- `additive_projection` 过高，说明虚频被 Lewis 酸运动主导。
- `reaction_projection` 过低，说明虚频未能描述目标反应坐标。
- 不通过时写入降级原因，并可强制 IRC 验证。

该逻辑对三种模式都仍然有效，应保留为 LA 模式下的 mandatory quality gate。

### 2.6 强度 registry

`rph_core/lewis_acid/acidity.py` 已提供 registry 框架：

```python
RegistryEntry(
    lewis_acid="AlCl3",
    acetone_affinity_kcal_mol=None,
    thf_affinity_kcal_mol=None,
    acetonitrile_affinity_kcal_mol=None,
    lumo_energy_ev=None,
    preferred_coordination_number=None,
    formal_charge=0,
    counterion_type="none",
)
```

当前不足：

- DFT 计算 runner 尚未实现。
- 没有 `acetone_delta_q_c`。
- 没有 charge method provenance。
- 没有 artifact 路径。
- 没有 small-molecule cache 联动。

因此目前可称为 Lewis 酸 registry 基础设施，而不是完整强度拟合系统。

## 3. 目标架构

### 3.1 双层设计

建议把 Lewis 酸体系明确拆成两层：

```text
RPH Lewis Acid System

Part A: Strength calibration
  位置: 小分子模块和 rph_core/lewis_acid
  输入: acetone + Lewis acid surrogate 或真实 Lewis acid
  输出: registry entry, calibration cache, delta_q_c

Part B: Pipeline additive mode
  位置: orchestrator, S1, S2, S3
  输入: none/strong/weak/custom mode
  输出: additive geometry, coordination trace, quality flags, pipeline metadata
```

Part A 回答“这个 Lewis 酸有多强”。  
Part B 回答“这个 Lewis 酸模式下的反应路径是否可信”。

### 3.2 数据流

```text
config/defaults.yaml
    |
    v
LewisAcidModeResolver
    |
    +--> LewisAcidAdditive
    |        |
    |        +--> S1 geometry SMILES append
    |        +--> S2 coordination trace
    |        +--> S3 imaginary mode projection
    |
    +--> LewisAcidCalibrationRunner
             |
             +--> cached free acetone reference charge
             +--> generated LA fragment geometry
             +--> acetone-LA complex
             +--> xTB pre-optimization and geometry filtering
             +--> charge extraction
             +--> config/lewis_acid_registry.json
             +--> SmallMolecules/lewis_acid_calibration/
```

## 4. 三模式配置方案

### 4.1 推荐配置结构

```yaml
lewis_acid:
  mode: none  # none | strong_single_site | weak_multi_site | custom_surrogate | legacy_licl

  modes:
    none:
      enabled: false
      surrogate: null
      strength_class: none

    strong_single_site:
      enabled: true
      surrogate: BF3
      coordination_model: single_site
      strength_class: strong
      s2_quality_gate:
        center_ligand_max: 2.6
        center_o_max: 2.3
        coordination_switch_allowed: false
      s3_quality_gate:
        reaction_mode_projection_min: 0.35
        additive_mode_projection_max: 0.30

    weak_multi_site:
      enabled: true
      surrogate: MgCl2
      coordination_model: multi_site
      strength_class: weak
      s2_quality_gate:
        center_ligand_max: 3.5
        center_o_max: 3.0
        coordination_switch_allowed: true
        top_n_coordination_sites: 2
      s3_quality_gate:
        reaction_mode_projection_min: 0.35
        additive_mode_projection_max: 0.20

    custom_surrogate:
      enabled: true
      surrogate: LiCl
      coordination_model: custom
      strength_class: custom

  surrogate_registry:
    BF3:
      name: BF3
      smiles: "B(F)(F)F"
      charge: 0
      multiplicity: 1
      center_element: B
      ligand_elements: [F]
      default_mode: strong_single_site

    MgCl2:
      name: MgCl2
      smiles: "Cl[Mg]Cl"
      charge: 0
      multiplicity: 1
      center_element: Mg
      ligand_elements: [Cl]
      default_mode: weak_multi_site
```

### 4.2 兼容旧配置

需要保留旧配置路径：

```yaml
lewis_acid:
  enabled: true
  surrogate: BF3
```

兼容解析规则：

| 旧配置 | 新解析结果 |
|---|---|
| `enabled: false` | `mode: none` |
| `enabled: true`, `surrogate: BF3` | `mode: strong_single_site` |
| `enabled: true`, `surrogate: MgCl2` | `mode: weak_multi_site` |
| `enabled: true`, 其他 surrogate | `mode: custom_surrogate` |

旧字段不应立即删除，只在日志中给出 debug 或 warning 级兼容提示。

### 4.3 模式解析器

建议新增：

```text
rph_core/lewis_acid/mode_resolver.py
```

核心 API：

```python
@dataclass
class LewisAcidModeSpec:
    mode: str
    enabled: bool
    surrogate_name: Optional[str]
    surrogate_smiles: Optional[str]
    charge: int
    multiplicity: int
    coordination_model: str
    strength_class: str
    center_element: Optional[str]
    ligand_elements: tuple[str, ...]
    s2_quality_gate: dict[str, Any]
    s3_quality_gate: dict[str, Any]
    calibration_required: bool = False


def resolve_lewis_acid_mode(config: dict[str, Any]) -> LewisAcidModeSpec:
    ...
```

解析器职责：

- 合并 root config、mode config、surrogate registry 和 legacy overrides。
- 校验 mode 是否已知。
- 校验 surrogate 是否存在。
- 校验 surrogate SMILES 是否可解析。
- 输出供 orchestrator 和 calibration runner 使用的稳定 dataclass。

失败策略：

- 配置错误 fail-fast。
- 缺少可选强度 registry 时 degrade gracefully。

## 5. Lewis 酸强度拟合方案

### 5.1 核心描述符

主描述符：

```text
acetone_delta_q_c = q_c_complex - q_c_free
```

字段定义：

| 字段 | 含义 |
|---|---|
| `acetone_q_c_free` | 游离丙酮羰基碳电荷 |
| `acetone_q_c_complex` | 丙酮-LA 配合物中羰基碳电荷 |
| `acetone_delta_q_c` | 配位导致的羰基碳电荷变化 |
| `charge_method_selected` | 实际使用的电荷方法 |
| `charge_methods_available` | 成功解析到的电荷方法列表 |

化学解释：

- `delta_q_c` 越正，羰基碳越缺电子。
- 对以羰基接受配位为标准探针的 Lewis 酸，`delta_q_c` 可作为电子吸引能力的统一尺度。
- 该指标是局部电子响应描述符，不直接等价于反应能垒或选择性。

### 5.2 电荷方法优先级

建议沿用配置式优先级：

```yaml
lewis_acid_calibration:
  charge_priority: ["NBO", "CM5", "MULLIKEN"]
```

解析策略：

1. 优先解析 NBO。
2. 如果没有 NBO，尝试 CM5。
3. 如果没有 CM5，尝试 Mulliken。
4. 如果全部缺失，写入 warning，`acetone_delta_q_c` 为 `null`，不让整个 pipeline 崩溃。

所有可用方法都应保存，避免后续复查时无法判断方法差异。

### 5.3 标准化强度分数

除原始 `acetone_delta_q_c` 外，可提供版本化归一分数：

```text
strength_score = clamp((delta_q_c - delta_q_none) / (delta_q_bf3 - delta_q_none), 0, 1)
```

初始锚点：

| 锚点 | 建议值 |
|---|---|
| `none` | `delta_q_none = 0.0`, `strength_score = 0.0` |
| `MgCl2` | 弱 Lewis 酸参考点 |
| `BF3` | 强 Lewis 酸参考点，默认上限 |

注意：

- `strength_score` 是拟合层结果，必须带 `fit_version`。
- `acetone_delta_q_c` 是原始 DFT 描述符，应长期保留。
- 后续真实 Lewis 酸库扩展时，可重新拟合 `strength_score`，但不应覆盖原始 `delta_q_c`。

### 5.4 辅助描述符

建议同时记录：

| 字段 | 用途 |
|---|---|
| `center_o_distance` | 配位几何质量 |
| `preferred_coordination_number` | 表征单配位或多配位倾向 |
| `lumo_energy_ev` | 可选电子结构描述符 |
| `coordination_model` | `single_site` 或 `multi_site` |
| `calibration_warnings` | 记录降级原因 |

本修正版不再把游离 Lewis 酸优化和 binding energy 作为主流程。原因是当前标定目标是羰基碳电荷响应，主公式只需要 `q_c(free acetone)` 和 `q_c(complex)`。游离 Lewis 酸能量只服务于 binding energy，而 binding energy 更容易受到溶剂、熵、聚集态和构象采样影响，会增加计算成本但不直接改善 `delta_q_c`。

如果后续确实需要 binding energy，可作为可选扩展任务单独开启：

```yaml
lewis_acid_calibration:
  optional_descriptors:
    binding_energy: false
```

## 6. 小分子校准流程

### 6.1 新增 runner

建议新增：

```text
rph_core/lewis_acid/calibration.py
```

核心类：

```python
class LewisAcidCalibrationRunner:
    def ensure_calibrated(
        self,
        lewis_acid_name: str,
        config: dict[str, Any],
        cache_root: Path,
    ) -> LewisAcidCalibrationRecord:
        ...
```

主要职责：

- 读取 mode spec 和 surrogate registry。
- 计算或复用游离丙酮的参考电荷；该步骤全局缓存，不随每个 Lewis 酸重复。
- 从 registry/template/RDKit 生成 Lewis 酸片段初始几何；不把游离 Lewis 酸 DFT 优化作为主流程必需步骤。
- 生成丙酮-LA 配合物初猜。
- 对候选配合物执行 xTB 预优化。
- 执行几何过滤。
- 选择最优有效配合物。
- 对最优配合物执行 Gaussian 单点和电荷分析。
- 提取配合物中羰基碳电荷。
- 计算 `delta_q_c = q_c(complex) - q_c(free)`。
- 写入 calibration cache。
- 更新 `config/lewis_acid_registry.json`。

所有 QC 调用必须通过 `rph_core.utils.qc_interface.py`，不能在 calibration runner 中直接调用 `subprocess.run`。

### 6.1.1 ZnCl2 标定修正版 8 步流程

针对 ZnCl2，推荐的最小有效标定流程如下：

1. 准备游离丙酮参考：如果 `SmallMolecules/lewis_acid_calibration/acetone/free/` 中已有结构和 `q_c(free)`，直接复用；只有缓存缺失时才计算。
2. 准备 ZnCl2 片段初猜：从内置模板或 registry geometry 读取线性或近线性 `Cl-Zn-Cl`，不把游离 ZnCl2 DFT 优化作为主流程步骤。
3. 构建 acetone-ZnCl2 配合物初猜：使用 3 个 O-Zn 初始距离乘以 12 个旋转，生成 36 个候选。
4. 对 36 个候选执行 xTB 预优化。
5. 做几何过滤：检查 C=O 拓扑、Zn-O 配位、Cl-Zn-Cl 骨架、原子重叠和明显解离。
6. 从通过过滤的候选中选择最低 xTB 能量结构作为 `best_complex.xyz`。
7. 仅对 `best_complex.xyz` 执行 Gaussian 单点和电荷提取，电荷优先级为 `NBO -> CM5 -> Mulliken`。
8. 计算 `delta_q_c = q_c(complex) - q_c(free)`，写入 `config/lewis_acid_registry.json` 的 v2 schema，并缓存到 `SmallMolecules/lewis_acid_calibration/ZnCl2/`。

该流程把计算量集中在真正影响 `delta_q_c` 的配合物构型选择和电荷提取上。游离丙酮只作为全局 reference；游离 ZnCl2 不进入主标定链路。

### 6.2 小分子缓存结构

建议缓存布局：

```text
SmallMolecules/
  lewis_acid_calibration/
    acetone/
      free/
        molecule_min.xyz
        finalDFT/
        charges.json
        cache_meta.json
    ZnCl2/
      acetone_complex/
        candidates/
          candidate_001_xtb.xyz
          candidate_002_xtb.xyz
          ...
        best_complex.xyz
        gaussian_sp/
        charges.json
        calibration.json
        cache_meta.json
```

缓存签名必须包含：

- Lewis 酸名称
- surrogate SMILES
- mode
- charge
- multiplicity
- opt method 和 basis
- sp method 和 basis
- solvent
- charge priority
- calibration schema version
- complex builder version

这样可以避免不同模式、不同理论级别或不同电荷方法之间错误复用缓存。

缓存根目录按真实或 surrogate Lewis 酸名称区分。例如 ZnCl2 标定写入：

```text
SmallMolecules/lewis_acid_calibration/ZnCl2/
```

### 6.3 丙酮-LA 配合物初猜

构建步骤：

1. 获取游离丙酮优化结构和 `q_c(free)`。若缓存已存在，直接复用。
2. 定位羰基 O 和羰基 C。
3. 从配置或模板获取 Lewis 酸片段几何，例如 ZnCl2 的线性或近线性 Cl-Zn-Cl 初始几何。
4. 沿羰基 O 的配位方向放置 Lewis 酸中心。
5. 生成多个 O-center 初始距离。
6. 围绕 C=O 轴旋转生成多个构象。
7. 对 ZnCl2 等二卤化物保留 Cl-Zn-Cl 骨架方向，并允许围绕 O-Zn 轴旋转。

建议配置：

```yaml
lewis_acid_calibration:
  complex_builder:
    rotations: 12
    max_candidates: 36
    center_o_distance_grid:
      BF3: [1.60, 1.75, 1.90]
      MgCl2: [2.00, 2.20, 2.40]
      ZnCl2: [2.00, 2.20, 2.40]
```

ZnCl2 标定的标准候选数：

```text
3 个 O-Zn 初始距离 x 12 个绕 O-Zn/C=O 相关轴旋转 = 36 个候选
```

有效性过滤：

- xTB 预优化收敛。
- acetone C=O 拓扑未改变。
- Zn 中心仍与羰基 O 配位。
- 中心-配体骨架未严重解离。
- 结构没有明显原子重叠。

选择策略：

- 优先选择通过几何过滤的最低 xTB 能量结构。
- 若存在多个近简并结构，记录 ensemble metadata。
- 对弱多配位模式，不强制单一 O 位点不变，但需要记录配位位点。

最优结构选出后，只对该结构执行 Gaussian 单点和电荷分析，不对全部 36 个候选执行 Gaussian。

### 6.4 推荐 QC 协议

配置应从主理论级别继承：

```yaml
lewis_acid_calibration:
  enabled: true
  reference_molecule: acetone
  charge_priority: ["NBO", "CM5", "MULLIKEN"]
  free_lewis_acid_optimization: false
  theory:
    xtb_preopt:
      inherit_from: theory.preoptimization
    sp:
      inherit_from: theory.single_point
    charge_analysis:
      enabled: true
```

不应硬编码：

- 方法
- 基组
- 溶剂
- 可执行文件路径
- charge method
- 资源数量
- 模板字符串

## 7. Registry v2 数据契约

### 7.1 目标 JSON

```json
{
  "schema_version": "la_registry_v2",
  "theory": {
    "opt_method": "B3LYP",
    "opt_basis": "def2-SVP",
    "sp_method": "wB97X-D3BJ",
    "sp_basis": "def2-TZVPP",
    "solvent": "acetone"
  },
  "calibration": {
    "reference_molecule": "acetone",
    "descriptor": "acetone_delta_q_c",
    "fit_version": "acetone_charge_response_v1",
    "charge_priority": ["NBO", "CM5", "MULLIKEN"]
  },
  "acids": {
    "ZnCl2": {
      "lewis_acid": "ZnCl2",
      "mode": "weak_multi_site",
      "surrogate_smiles": "Cl[Zn]Cl",
      "acetone_q_c_free": 0.6123,
      "acetone_q_c_complex": 0.6815,
      "acetone_delta_q_c": 0.0692,
      "charge_method_selected": "NBO",
      "charge_methods_available": ["NBO", "MULLIKEN"],
      "strength_score": 0.46,
      "strength_class": "weak",
      "center_o_distance": 2.18,
      "preferred_coordination_number": 2,
      "formal_charge": 0,
      "counterion_type": "chloride",
      "artifacts": {
        "free_acetone": "SmallMolecules/lewis_acid_calibration/acetone/free",
        "acetone_complex": "SmallMolecules/lewis_acid_calibration/ZnCl2/acetone_complex"
      },
      "warnings": ["free_lewis_acid_optimization_skipped"]
    }
  }
}
```

### 7.2 RegistryManager 行为

`RegistryManager.load()` 应支持：

- v1 registry 正常读取。
- v2 registry 正常读取。
- 缺失的新字段填为 `None`。
- corrupt entry 跳过并记录 warning。
- 缺失 registry 返回空 dict，不中断主流程。

`RegistryManager.save()` 应支持：

- 默认写出 `la_registry_v2`。
- 保留已有 v1 字段。
- 写入 `theory` 与 `calibration` 元数据。
- 写入 artifact 相对路径。

## 8. Pipeline 集成方案

### 8.1 S0 条件检测

当前 `detect_lewis_acid()` 只返回 bool。建议升级为结构化结果：

```python
@dataclass
class LewisAcidDetectionResult:
    use_la: bool
    canonical_name: Optional[str]
    aliases_matched: tuple[str, ...]
    suggested_mode: str
    suggested_surrogate: Optional[str]
    confidence: str
```

映射应配置化：

```yaml
lewis_acid_detection:
  mapping:
    BF3:
      mode: strong_single_site
      surrogate: BF3
    BCl3:
      mode: strong_single_site
      surrogate: BF3
    AlCl3:
      mode: strong_single_site
      surrogate: BF3
    MgCl2:
      mode: weak_multi_site
      surrogate: MgCl2
    ZnCl2:
      mode: weak_multi_site
      surrogate: MgCl2
```

当真实 Lewis 酸尚未校准时，可先映射到强/弱 surrogate；校准完成后再使用 registry 中的真实 `strength_score`。

### 8.2 S1

保留现有原则：

- canonical product SMILES 不变。
- geometry product SMILES 拼接 surrogate。
- additive 原子位于 XYZ 末尾。
- S1 sidecar 写入原子映射与 additive metadata。

新增 sidecar 字段：

```json
{
  "lewis_acid": {
    "enabled": true,
    "mode": "strong_single_site",
    "surrogate_name": "BF3",
    "surrogate_smiles": "B(F)(F)F",
    "coordination_model": "single_site",
    "strength_class": "strong",
    "strength_score": 1.0,
    "acetone_delta_q_c": 0.1218,
    "additive_elements": ["B", "F", "F", "F"],
    "additive_atom_indices": [40, 41, 42, 43],
    "center_atom_index": 40,
    "organic_atom_count": 40
  }
}
```

### 8.3 S2

S2 trace 应从 Li 命名升级到 center 命名。

兼容映射：

| 旧字段 | 新字段 | 策略 |
|---|---|---|
| `li_cl_distances` | `center_ligand_distances` | 新旧都写 |
| `li_o_distances` | `center_o_distances` | 新旧都写 |
| `li_coordination_sites` | `coordination_sites` | 新旧都写 |
| `metal_cl_max` | `center_ligand_max` | 旧字段作为 alias |
| `metal_o_max` | `center_o_max` | 旧字段作为 alias |

S2 trace v2 示例：

```json
{
  "schema": "la_coordination_trace_v2",
  "mode": "weak_multi_site",
  "center_atom_index": 42,
  "center_element": "Mg",
  "ligand_atom_indices": [41, 43],
  "frames": [
    {
      "step": 0,
      "center_ligand_min": 2.21,
      "center_o_min": 2.08,
      "coordination_sites": [12, 19]
    }
  ],
  "coordination_switch_count": 1,
  "coordination_switch_allowed": true
}
```

### 8.4 S3

`la_mode_projection.json` 建议增加：

```json
{
  "schema": "la_mode_projection_v2",
  "la_enabled": true,
  "mode": "strong_single_site",
  "surrogate_name": "BF3",
  "center_atom_index": 40,
  "coordination_model": "single_site",
  "strength_class": "strong",
  "mode_index": 0,
  "frequency_cm-1": -432.1,
  "additive_projection": 0.08,
  "reaction_projection": 0.62,
  "passed": true,
  "reason": "",
  "irc_forced": false
}
```

不建议让 S3 LA projection 失败直接 crash。更稳妥的策略：

- 标记 `DEGRADED`。
- 写入 `LewisAcidQualityFlags`。
- 根据 `require_irc_for_flagged_ts` 强制 IRC 或提高验证等级。

### 8.5 pipeline.state 元数据

建议写入 compact metadata：

```json
{
  "lewis_acid": {
    "enabled": true,
    "mode": "strong_single_site",
    "surrogate_name": "BF3",
    "surrogate_smiles": "B(F)(F)F",
    "strength_class": "strong",
    "strength_score": 1.0,
    "acetone_delta_q_c": 0.1218,
    "additive_atom_indices": [40, 41, 42, 43],
    "center_atom_index": 40,
    "organic_atom_count": 40,
    "quality_flags": {
      "s2_passed": true,
      "s3_passed": true,
      "additive_projection": 0.08,
      "reaction_projection": 0.62
    }
  }
}
```

这些数据只作为外部 S4 和 ML 的输入元数据，不在 RPH 内部生成 S4 特征。

## 9. 文件级实施计划

| 文件 | 动作 | 目的 |
|---|---|---|
| `rph_core/lewis_acid/models.py` | 修改 | 增加 mode/calibration dataclass，添加通用 trace 命名 |
| `rph_core/lewis_acid/mode_resolver.py` | 新增 | 统一解析 `none/strong/weak/custom` |
| `rph_core/lewis_acid/acidity.py` | 修改 | registry v2，兼容 v1 |
| `rph_core/lewis_acid/calibration.py` | 新增 | 丙酮-LA 强度拟合 runner |
| `rph_core/utils/small_molecule_cache.py` | 修改 | 增加 LA calibration cache namespace |
| `rph_core/utils/small_molecule_catalog.py` | 修改 | 确保 acetone 作为标准探针可复用 |
| `rph_core/utils/charge_reader.py` | 修改 | 增加优先级电荷提取 API |
| `rph_core/utils/lewis_acid_detector.py` | 修改 | 返回结构化检测结果 |
| `rph_core/orchestrator.py` | 修改 | 使用 mode spec，写入更完整 LA metadata |
| `rph_core/steps/step2_retro/retro_scanner.py` | 修改 | 泛化 S2 trace 与 gate 命名 |
| `rph_core/steps/step3_opt/ts_optimizer.py` | 修改 | projection report 增加 mode metadata |
| `config/defaults.yaml` | 修改 | 增加三模式配置与 calibration 配置 |
| `tests/test_lewis_acid_modes.py` | 新增 | mode 解析测试 |
| `tests/test_lewis_acid_calibration.py` | 新增 | mocked calibration、registry 和 cache 测试 |

## 10. 分阶段路线

### Phase 1: 模式层稳定化

目标：

- 用户可以直接配置 `none`、`strong_single_site`、`weak_multi_site`。
- 旧 `enabled/surrogate` 配置继续可用。

任务：

1. 新增 `mode_resolver.py`。
2. 在 orchestrator 中用 `LewisAcidModeSpec` 替代散落的 raw config 解析。
3. 给 `LewisAcidAdditive` 增加 mode、strength_class、coordination_model 等字段。
4. 补 `test_lewis_acid_modes.py`。

验收：

- `none` 不拼接 surrogate。
- `strong_single_site` 拼接 BF3。
- `weak_multi_site` 拼接 MgCl2。
- 旧 BF3/MgCl2 配置能自动映射到新模式。

### Phase 2: S2/S3 泛化

目标：

- 消除 LiCl 专属语义对 BF3/MgCl2 的干扰。
- 保留旧字段兼容性。

任务：

1. 新增 `MetalCoordinationTrace` 或 `CenterCoordinationTrace`。
2. 保留 `LiCoordinationTrace` 作为 alias 或兼容 dataclass。
3. S2 同时写新旧 trace 字段。
4. S2 gate 使用 `center_ligand_max`、`center_o_max`。
5. S3 projection report 写入 mode metadata。

验收：

- BF3 走 B-F/B-O 检查。
- MgCl2 允许 coordination switch。
- 旧测试 `test_s2_la_trace.py` 仍通过。

### Phase 3: 小分子强度校准

目标：

- 在小分子模块中完成 acetone-LA 强度拟合闭环；主流程只计算 `q_c(free acetone)` 和 `q_c(complex)`，不再强制优化游离 Lewis 酸。

任务：

1. 新增 `LewisAcidCalibrationRunner`。
2. 复用或扩展 `SmallMoleculeCache`。
3. 建立全局 free acetone reference cache。
4. 为 ZnCl2 构建 36 个 acetone-ZnCl2 complex 候选。
5. 对候选执行 xTB 预优化、几何过滤和最优结构选择。
6. 对最优 complex 执行 Gaussian 单点和电荷提取。
7. 计算 `delta_q_c = q_c(complex) - q_c(free)`。
8. 写入 `calibration.json` 和 registry v2。

验收：

- mocked QC 下能生成 `acetone_delta_q_c`。
- 缓存命中时不重复计算。
- 缺失 NBO 时能 fallback。
- 全部 charge 缺失时记录 warning，不 crash。
- registry 中 `ZnCl2` 条目不要求 `free_lewis_acid` artifact。

### Phase 4: 条件自动映射

目标：

- 根据 cleaner_data 中 catalyst/additive 自动建议 mode。

任务：

1. 升级 `lewis_acid_detector.py`。
2. 增加 config-driven mapping。
3. 在 S0 summary 写入 canonical LA、suggested mode、confidence。

验收：

- BF3 条件自动映射 strong。
- MgCl2 条件自动映射 weak。
- 未识别条件保持 none 或用户显式配置。

### Phase 5: 文档与外部契约

目标：

- 明确 RPH 输出哪些 LA metadata。
- 明确 RPH_Postprocess 可消费哪些字段。

任务：

1. 更新 README.zh-CN.md 的 LA 配置示例。
2. 更新 `LEWIS_ACID_MODULE_DESIGN.md` 或添加链接到本文档。
3. 给 `pipeline.state`、`scan_profile.json`、`la_mode_projection.json` 写 schema 说明。

验收：

- 用户能按文档跑三种模式。
- 下游能稳定读取 LA 元数据。

## 11. 测试计划

### 11.1 单元测试

模式解析：

- `mode=none` 禁用 LA。
- `mode=strong_single_site` 解析为 BF3。
- `mode=weak_multi_site` 解析为 MgCl2。
- `mode=custom_surrogate` 需要显式 surrogate。
- 未知 mode 报清晰配置错误。

registry：

- v1 registry 可读取。
- v2 registry 可读取。
- 缺字段填 `None`。
- corrupt entry 记录 warning 并跳过。

calibration：

- free acetone 只计算一次。
- BF3 complex 生成 `delta_q_c`。
- MgCl2 complex 记录多配位 metadata。
- NBO 缺失时 fallback CM5 或 Mulliken。
- 所有 charge 缺失时 warning，不 crash。

S2：

- BF3 检查 center-ligand 和 center-O。
- MgCl2 允许配位切换。
- 旧 `li_o_max`、`li_cl_max` alias 仍工作。

S3：

- projection report 包含 mode。
- additive-dominated TS 设置 degraded flags。
- 缺 displacement data 不阻断主流程。

### 11.2 集成测试

| 场景 | 预期 |
|---|---|
| no LA | geometry SMILES 不拼 surrogate，不写 LA trace |
| BF3 strong | geometry 包含 BF3，S2/S3 LA metadata 存在 |
| MgCl2 weak | geometry 包含 MgCl2，配位切换可被允许 |
| auto-detected BF3 | S0 建议 strong mode |
| auto-detected ZnCl2 | 未校准时映射到 weak fallback |

### 11.3 提交前命令

```bash
python scripts/ci/check_imports.py rph_core
pytest tests/test_lewis_acid_modes.py -v
pytest tests/test_lewis_acid_calibration.py -v
pytest tests/test_s2_la_trace.py -v
pytest tests/test_s3_la_validator.py -v
```

注意：当前本地 Python 环境若缺少 RDKit，Lewis acid 相关测试可能在 collection 阶段因 `ModuleNotFoundError: No module named 'rdkit'` 失败。正式验收环境需要安装 RDKit，或将纯 dataclass 测试从 package-level orchestrator import 中解耦。

## 12. 迁移策略

### 12.1 保持兼容

必须继续支持：

- `lewis_acid.enabled`
- `lewis_acid.surrogate`
- `surrogate_overrides`
- `LewisAcidAdditive.metal_atom_index`
- `LiCoordinationTrace`
- `li_o_distances`
- `li_cl_distances`

### 12.2 渐进替换

第一阶段：

- 新增新字段。
- 同时写新旧字段。
- 旧配置不报警或仅 debug 提示。

第二阶段：

- 文档全部改用新模式。
- 测试覆盖新旧路径。

第三阶段：

- 对历史内部 mode 字符串给 warning。
- 等下游 RPH_Postprocess 确认兼容后，再考虑删除旧字段。

## 13. 风险与缓解

| 风险 | 影响 | 缓解 |
|---|---|---|
| BF3 配合物优化时解离 | 强度 descriptor 无效 | 加几何过滤和 candidate rejection |
| MgCl2 多配位浅势能面导致结果噪声 | 弱酸强度不稳定 | 多初猜采样并记录 ensemble metadata |
| 不同 charge 方法排序不一致 | 强度排序不稳 | 保存所有方法，`strength_score` 标注 fit_version |
| 不同理论级别复用 cache | descriptor 不一致 | cache signature 纳入 theory 和 charge priority |
| Li 命名残留 | 后续维护混乱 | 新增 center 命名并保留旧 alias |
| 自动检测过粗 | 真实 Lewis 酸映射错误 | config-driven mapping 加 confidence |
| 校准拖慢主流程 | 批量计算成本高 | calibration 默认缓存化，可预计算 |

## 14. 完成定义

该方案完成的判断标准：

1. 用户可通过 `none`、`strong_single_site`、`weak_multi_site` 配置三种 Lewis 酸模式。
2. 强酸模式默认 BF3，弱酸模式默认 MgCl2，且都能正确进入 S1/S2/S3。
3. 小分子模块能计算或复用 acetone-LA calibration。
4. registry v2 能保存 `acetone_delta_q_c`、`strength_score`、charge method 和 artifact。
5. S2 能写通用 center coordination trace。
6. S3 能写带 mode metadata 的 projection report。
7. 旧 LiCl 配置和旧 trace 字段保持兼容。
8. 在含 RDKit 的环境中，LA 相关单元测试和 import gate 通过。
