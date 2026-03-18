# S2 阶段 xTB Path Finder 集成改造方案

**文档版本**: v1.0  
**编制日期**: 2026-03-15  
**状态**: 待评审  

---

## 一、改造背景与目的

### 1.1 问题陈述

当前 S2 阶段使用 **xTB relaxed scan** 方法生成过渡态初猜，该方法通过在反应坐标上拉伸键长来寻找能量最高点。然而：

1. **拉伸生成方案效果差**：此前的 kinematic stretch + 中间体优化方案（forward_scan）已被用户反馈效果不佳
2. **当前 retro_scan 仍不理想**：改为直接 outward scan 后，用户仍反馈效果较差
3. **缺乏多样性**：目前只有单一扫描策略，缺乏更robust的路径搜索能力

### 1.2 改造目标

引入 **xTB 官方 Path Finder（元动力学反应路径搜索）** 作为 S2 的新主流程，以实现：

- 更可靠的 TS 初猜生成（基于 start → end 端点的路径搜索）
- 与 Q-Chem/Sella 的-growing string method (GSM) 类似的第一性原理方法
- 保留现有架构的最大复用，避免大规模重构

---

## 二、现状分析

### 2.1 当前 S2 执行流程

```
S1 (Product) → S2 (TS Guess) → S3 (TS Optimization) → S4 (Feature Extraction)
                        ↓
              RetroScanner.run_retro_scan()
                        ↓
              XTBInterface.scan() / XTBRunner.run_scan()
                        ↓
              产出: ts_guess.xyz, reactant_complex.xyz, dipolar_intermediate.xyz
```

### 2.2 当前代码架构

| 组件 | 位置 | 职责 |
|------|------|------|
| `RetroScanner` | `rph_core/steps/step2_retro/retro_scanner.py` | S2 协调器 |
| `XTBInterface` | `rph_core/utils/qc_interface.py` | xTB 封装 |
| `XTBRunner` | `rph_core/utils/xtb_runner.py` | xTB 执行器 |
| `ScanResult` | `rph_core/utils/data_types.py` | 扫描结果数据模型 |

### 2.3 下游依赖分析

| 消费者 | 依赖文件 | 依赖方式 |
|--------|----------|----------|
| S3 (Berny) | `ts_guess.xyz` | TS 优化起点 |
| S3 (QST2) | `reactant_complex.xyz` + `product.xyz` | QST2 端点 |
| Checkpoint | `ts_guess_xyz`, `reactant_xyz` | 恢复依据 |
| S4 (Mech Packager) | `reactant_complex.xyz` | 机理分析回退源 |
| S4 (Fragment) | `reactant_xyz` | 片段分割必需 |

**关键发现**：`reactant_complex.xyz` 在多处被当作"真实反应物"使用，而非仅仅是"扫描起点"。

---

## 三、xTB Path Finder 技术细节

### 3.1 官方 CLI 用法

```bash
xtb start.xyz --path end.xyz --input path.inp
```

**参数说明**：
- `start.xyz`：反应物/中间体结构
- `end.xyz`：产物结构
- `path.inp`：路径搜索参数

### 3.2 path.inp 推荐参数

```xcontrol
$path
   nrun=1          # 路径优化轮数
   npoint=25       # 初始路径点数
   anopt=10       # 精细优化步数
   kpush=0.003    # 推力（从反应物出发）
   kpull=-0.015   # 拉力（向产物方向）
   ppull=0.05     # 优化器拉力
   alp=1.2        # 高斯宽度
$end
```

### 3.3 输出文件

| 文件 | 含义 |
|------|------|
| `xtbpath_ts.xyz` | **预估过渡态**（关键输出） |
| `xtbpath_0.xyz` | 初始路径 |
| `xtbpath_*.xyz` | 各轮优化路径 |

### 3.4 输出解析

```
path  2 taken with   23 points.
estimated TS on file xtbpath_ts.xyz

forward  barrier (kcal) :    12.420
backward barrier (kcal) :    37.497
reaction energy  (kcal) :   -25.076
norm(g) at est. TS, point: 0.01615
```

---

## 四、架构设计

### 4.1 设计原则

1. **渐进式迁移**：不破坏现有契约，新增角色命名
2. **最大复用**：复用工件目录、sandbox、日志、错误处理
3. **可验证性**：每一步产出可追溯

### 4.2 新增数据模型

```python
# rph_core/utils/data_types.py

@dataclass
class PathSearchResult:
    """xTB 路径搜索结果"""
    success: bool = False
    path_xyz_files: List[Path] = None  # 完整路径轨迹
    ts_guess_xyz: Optional[Path] = None  # 预估 TS (xtbpath_ts.xyz)
    path_log: Optional[Path] = None  # 日志文件
    
    # 能量信息
    barrier_forward_kcal: Optional[float] = None
    barrier_backward_kcal: Optional[float] = None
    reaction_energy_kcal: Optional[float] = None
    
    # TS 质量指标
    estimated_ts_point: Optional[int] = None
    gradient_norm_at_ts: Optional[float] = None
    
    error_message: Optional[str] = None
```

### 4.3 S2 → S3 契约演进

| 阶段 | ts_guess.xyz | dipolar_intermediate.xyz | reactant_complex.xyz |
|------|---------------|-------------------------|---------------------|
| 现有 | 扫描能量最高点 | 扫描起点 | 扫描终点（同名） |
| 过渡 | xtbpath_ts.xyz | path 起点 | dipolar_intermediate 的兼容别名 |
| 目标 | xtbpath_ts.xyz | path 起点 | 退役 |

**元数据声明**：
```json
{
  "generation_method": "xtb_path",
  "start_structure_role": "dipolar_intermediate",
  "reactant_complex_is_alias": true
}
```

### 4.4 S3 Artifact Resolver

在 Step3 入口增加轻量解析层：

```python
def resolve_s3_inputs(s2_dir: Path, s3_config: dict) -> S3Inputs:
    """
    解析 S3 输入，优先使用新角色命名，回退到 legacy
    """
    # 1. 检查元数据声明
    metadata = load_s2_metadata(s2_dir)
    
    # 2. 优先使用新角色
    if metadata.get("start_structure_role") == "dipolar_intermediate":
        return S3Inputs(
            ts_guess=ts_guess_xyz,
            reactant=dipolar_intermediate_xyz,
            product=product_xyz,
            source="dipolar_intermediate"
        )
    
    # 3. 回退到 legacy
    return S3Inputs(
        ts_guess=ts_guess_xyz,
        reactant=reactant_complex_xyz,
        product=product_xyz,
        source="legacy_reactant_complex"
    )
```

---

## 五、详细实施计划

### 5.1 第一阶段：基础设施（预计 2h）

#### 5.1.1 新增 PathSearchResult
- **文件**: `rph_core/utils/data_types.py`
- **内容**: 添加 `PathSearchResult` dataclass

#### 5.1.2 新增 XTBRunner.path()
- **文件**: `rph_core/utils/xtb_runner.py`
- **新增方法**: `def run_path(...)`
- **复用**: 
  - `_verify_executable()` - 可执行文件查找
  - `_run_command()` - 命令执行
  - 工作目录/sandbox 处理
  - 日志与错误处理

```python
def run_path(
    self,
    start_xyz: Path,
    end_xyz: Path,
    workdir: Path,
    nrun: int = 1,
    npoint: int = 25,
    anopt: int = 10,
    kpush: float = 0.003,
    kpull: float = -0.015,
    ppull: float = 0.05,
    alp: float = 1.2,
    gfn_level: int = 2,
    solvent: Optional[str] = None,
    charge: int = 0,
    uhf: int = 0,
    etemp: Optional[float] = None,
) -> PathSearchResult:
    """运行 xTB 路径搜索"""
    # 1. 写 path.inp
    # 2. 执行 xtb start.xyz --path end.xyz --input path.inp
    # 3. 解析输出，提取 xtbpath_ts.xyz 和能量信息
    # 4. 返回 PathSearchResult
```

#### 5.1.3 新增 XTBInterface.path()
- **文件**: `rph_core/utils/qc_interface.py`
- **新增方法**: `def path(...)` - 对外统一接口

### 5.2 第二阶段：S2 改造（预计 2h）

#### 5.2.1 新增 S2 配置块
- **文件**: `config/defaults.yaml`
- **内容**:
```yaml
step2:
  mode: path_search  # 或继续用 scan
  path_search:
    enabled: true
    nrun: 1
    npoint: 25
    anopt: 10
    kpush: 0.003
    kpull: -0.015
    ppull: 0.05
    alp: 1.2
    write_legacy_reactant_alias: true  # 兼容模式
```

#### 5.2.2 修改 RetroScanner
- **文件**: `rph_core/steps/step2_retro/retro_scanner.py`
- **修改**: 在 `run_retro_scan()` 中：
  1. 检测 `path_search.enabled`
  2. 若启用：调用 `XTBInterface.path()`
  3. 产出 `xtbpath_ts.xyz` → `ts_guess.xyz`
  4. 保留 `dipolar_intermediate.xyz` 作为起点
  5. 写 `reactant_complex.xyz` 作为兼容别名

#### 5.2.3 修改 scan_profile.json
- 升级为 `reaction_path_profile.json`（或兼容保留）
- 字段变更：
```json
{
  "generation_method": "xtb_path",
  "path_parameters": {
    "nrun": 1,
    "npoint": 25,
    ...
  },
  "energies": {
    "barrier_forward_kcal": 12.4,
    "barrier_backward_kcal": 37.5,
    "reaction_energy_kcal": -25.1
  }
}
```

### 5.3 第三阶段：S3 适配（预计 1h）

#### 5.3.1 新增 Artifact Resolver
- **文件**: `rph_core/steps/step3_opt/input_resolver.py`（新增）
- **职责**:
  - 读取 S2 元数据
  - 解析 start_structure 角色
  - 返回正确的输入路径

#### 5.3.2 修改 run_step3 入口
- **文件**: `rph_core/steps/runners.py`
- **修改**: 在调用 S3 前先经过 Resolver

### 5.4 第四阶段：Checkpoint 与 Orchestrator（预计 1h）

#### 5.4.1 修改 Checkpoint Manager
- **文件**: `rph_core/utils/checkpoint_manager.py`
- **修改**:
  - 允许从 `dipolar_intermediate_xyz` 恢复 `reactant_xyz`
  - 记录 `generation_method=xtb_path`

#### 5.4.2 修改 Orchestrator
- **文件**: `rph_core/orchestrator.py`
- **修改**: 
  - 日志文案改为 role-aware
  - Resume 逻辑支持新模式

### 5.5 第五阶段：测试（预计 2h）

#### 5.5.1 新增单元测试
- `tests/test_xtb_path.py` - XTBRunner.path() 测试
- `tests/test_path_resolver.py` - Artifact Resolver 测试

#### 5.5.2 修改现有测试
- 更新所有硬编码 `reactant_complex.xyz` 的测试
- 添加兼容模式验证

---

## 六、风险评估

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|----------|
| xTB --path 执行失败 | 中 | 高 | 添加 fallback 到原有 scan |
| 下游误用 alias | 高 | 高 | 元数据明确标记 |
| Checkpoint 不兼容 | 中 | 中 | 保留双重恢复逻辑 |
| 测试大规模失败 | 中 | 中 | 分阶段验证 |

---

## 七、验收标准

1. **功能验收**：
   - `xtb start.xyz --path end.xyz` 能正确产出 `xtbpath_ts.xyz`
   - S2 产出文件满足 S3 入口要求
   
2. **回归验收**：
   - 现有 S3/QST2/Checkpoint/S4 功能不受影响
   - 现有测试通过率 ≥ 95%

3. **文档验收**：
   - 更新相关 AGENTS.md
   - 更新 config/defaults.yaml 中文注释

---

## 八、附录

### A. 相关文件清单

| 类别 | 文件 | 操作 |
|------|------|------|
| 数据模型 | `rph_core/utils/data_types.py` | 修改 |
| 执行层 | `rph_core/utils/xtb_runner.py` | 修改 |
| 封装层 | `rph_core/utils/qc_interface.py` | 修改 |
| S2 逻辑 | `rph_core/steps/step2_retro/retro_scanner.py` | 修改 |
| 配置 | `config/defaults.yaml` | 修改 |
| S3 入口 | `rph_core/steps/runners.py` | 修改 |
| S3 解析器 | `rph_core/steps/step3_opt/input_resolver.py` | 新增 |
| Checkpoint | `rph_core/utils/checkpoint_manager.py` | 修改 |
| Orchestrator | `rph_core/orchestrator.py` | 修改 |
| 测试 | `tests/test_xtb_path.py` | 新增 |
| 测试 | `tests/test_path_resolver.py` | 新增 |

### B. 命令行验证

```bash
# 本地验证 xTB --path 可用性
xtb start.xyz --path end.xyz --input path.inp

# 检查输出
ls -la xtbpath*.xyz
```

### C. 参考资料

- [xTB Reaction Path Methods Documentation](https://xtb-docs.readthedocs.io/en/latest/path.html)
- [xcontrol $path section](https://github.com/grimme-lab/xtb/blob/main/man/xcontrol.7.adoc)
- [Growing String Method (GSM)](https://xtb-docs.readthedocs.io/en/latest/gsm.html)

---

**文档结束**
