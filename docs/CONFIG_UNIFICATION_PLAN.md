# RPH 配置层统一修改计划（v2.0）

> **版本**: v2.0 | **日期**: 2026-06-30 | **状态**: 待审核

---

## 一、问题汇总

### 1.1 配置碎片化

同一台机器的同一个资源参数分散在 **7 个 YAML 区段 + 15 处代码默认值** 中：

| 区段（YAML key 路径） | 用途 | 当前值 | 全局资源 |
|----------------------|------|--------|---------|
| `resources.nproc/mem` | 全局基数 | 16 / 32GB | **单源真值**（目标） |
| `theory.optimization.nproc/mem` | Gaussian OPT | 16 / (代码 16GB) | 应有继承 |
| `theory.single_point.nproc` | ORCA SP | 16 | 应有继承 |
| `step1.crest.threads` | CREST 构象搜索 | 20 | 应有继承 |
| `step1.conformer_search.common.threads` | CREST 第二回退 | 20 | 应有继承 |
| `step2.scan.nproc` | xTB constrained scan | 16 | 应有继承 |
| `step2.path_search.nproc` | xTB path search | 16 | 应有继承 |
| `step2.xtb_settings.nproc` | xTB GFN2 scan | 16 | 应有继承 |
| `intra_reaction_parallel.total_cores` | 并行预算池 | 20 | 应有自动检测 |
| `intra_reaction_parallel.total_mem` | 并行预算池 | 60GB | 应有自动检测 |
| `intra_reaction_parallel.s1.crest_cores` | S1 CREST 单 job | 20 | 应预算分配 |
| `intra_reaction_parallel.s1.opt_cores_per_job` | S1 OPT 单 job | 20 | 应预算分配 |
| `intra_reaction_parallel.s1.sp_cores_per_job` | S1 SP 单 job | 20 | 应预算分配 |
| `intra_reaction_parallel.s3.ts_opt_cores` | S3 TS 单 job | 20 | 应预算分配 |
| `intra_reaction_parallel.s3.sp_cores_per_job` | S3 SP 单 job | 20 | 应预算分配 |

### 1.2 代码默认值冲突

| 严重程度 | 文件 | Key | 代码默认 | 全局资源 | 影响 |
|---------|------|-----|---------|---------|------|
| 🔴 BUG | `gau_xtb_interface.py` | `%nprocshared` | **永远写 1** | 任意 | Gau_XTB 永远单核 |
| 🔴 BUG | `runners.py` | GauXTBInterface nproc | **硬编码 1** | - | 同上 |
| 🟡 | `engine.py` | `theory.optimization.mem` | **"16GB"** | 32GB | S1 DFT 内存砍半 |
| 🟡 | `qc_interface.py` | `run_gaussian_calc mem` | **"48GB"** | 32GB | 独立路径虚高 50% |
| 🟡 | `orca_interface.py` | maxcore safety | **0.8** | 0.65 | ORCA 内存高 23% |
| 🟡 | `resource_utils.py` | maxcore safety | **0.8** | 0.65 | |
| 🟡 | `ts_optimizer.py` | SP xtb nproc | **8** | 16 | 核数砍半 |
| 🟡 | `xtb_runner.py` | resources.nproc 回退 | **1** | 16 | xTB 默认单核 |
| 🟡 | `qc_interface.py:XTBInterface` | nproc 构造器默认 | **1** | 16 | |
| 🟡 | `qc_interface.py:CRESTInterface` | nproc 构造器默认 | **1** | 16 | |
| 🟡 | `intra_reaction_scheduler.py` | 多个并行 lane 默认 | 4/8/12/"24GB" | 16/32GB | 严重不符 |

### 1.3 可执行文件路径问题

- 仓库中共 **6 个活跃 YAML 配置文件**，全部硬编码 `/home/xieningke/...` 路径（已手动修复到 `/opt/software/...`，但无自动发现机制）
- 6 种软件的发现逻辑分散在不同模块，行为不一致
- Shermo 和 Multiwfn **完全没有回退机制**
- 无开机自动发现机制

**受影响配置文件**：`defaults.yaml`, `defaults_alcl3_test.yaml`, `defaults_la_test.yaml`, `defaults_la_test_dataset.yaml`, `defaults_la_test_fast.yaml`, `validation.yaml`

---

## 二、实施顺序（核心）

按 **先修 BUG → 再统一配置 → 最后加自动发现** 的顺序分阶段执行，每阶段独立可测试。

---

### 阶段 1: 修复致命 BUG

#### 1.1 Gau_XTB 单核 BUG

**问题**：`gau_xtb_interface.py` 的 `write_input_file()` 始终写 `%nprocshared=1`，无视构造器 `nproc` 参数。

```python
# ❌ 当前 (gau_xtb_interface.py: write_input_file)
f.write("%nprocshared=1\n")     # ← 永远单核

# ✅ 修复
f.write(f"%nprocshared={self.nproc}\n")
```

**连带修复**：`runners.py` 的调用点 `GauXTBInterface(config=..., nproc=1)` → 改为从 `resources.nproc` 获取。

**涉及文件**：`gau_xtb_interface.py`, `runners.py`

#### 1.2 统一 maxcore safety 系数

**问题**：YAML 配 `orca_maxcore_safety: 0.65`，但代码默认 `0.8`，导致 ORCA maxcore 高 23%。

**修改**：`orca_interface.py` 和 `resource_utils.py` 的代码默认值从 `0.8` 改为 `0.65`。

```python
# resource_utils.py
def calc_orca_maxcore(mem, nproc, safety_factor=0.65):  # 0.8→0.65

# orca_interface.py
safety_factor = res_cfg.get('orca_maxcore_safety', 0.65)  # 0.8→0.65
```

**涉及文件**：`resource_utils.py`, `orca_interface.py`

#### 1.3 统一 theory.optimization.mem 默认值

**问题**：`engine.py` 从 `theory.optimization.mem` 读 mem 时默认值为 `"16GB"`，是 YAML 的 32GB 和 ts_optimizer.py 默认值（"32GB"）的一半。

```python
# engine.py
mem=self.theory_opt.get('mem', '32GB')  # '16GB'→'32GB'
```

**涉及文件**：`engine.py`（四处）

#### 1.4 统一 run_gaussian_calculation mem 默认值

**问题**：`qc_interface.py` 中 `run_gaussian_calculation()` 默认 `"48GB"`，其他路径都是 32GB。

```python
mem = config.get('mem', '32GB')  # '48GB'→'32GB'
```

**涉及文件**：`qc_interface.py`

---

### 阶段 2: 统一配置继承层级

#### 2.1 Override 语义定义

| YAML 值 | 语义 | 示例 |
|---------|------|------|
| `null`（或不存在） | 继承自 `resources.*` + budget 分配 | `step1.crest.threads: null` |
| 数值/字符串 | **显式 override**，不继承不计算 | `step1.crest.threads: 12` |

**所有区段的 nproc/mem/threads 显式数值都是 override**，`null` 才是继承。

#### 2.2 清理 `defaults.yaml`（YAML key 路径）

将下列区段的 nproc/mem/threads **全部改为 `null`**：

| YAML key 路径 | 当前值 | 改为 |
|---------------|--------|------|
| `theory.optimization.nproc` | 16 | `null` |
| `theory.optimization.mem` | 32GB | `null` |
| `theory.single_point.nproc` | 16 | `null` |
| `step1.crest.threads` | 20 | `null` |
| `step1.conformer_search.common.threads` | 20 | `null` |
| `step2.scan.nproc` | 20 | `null` |
| `step2.path_search.nproc` | 20 | `null` |
| `step2.xtb_settings.nproc` | 20 | `null` |
| `step2.crest_rescue.threads` | 20 | `null` |
| `intra_reaction_parallel.total_cores` | 20 | `null` |
| `intra_reaction_parallel.total_mem` | "60GB" | `null` |
| `intra_reaction_parallel.s1.opt_cores_per_job` | 20 | `null` |
| `intra_reaction_parallel.s1.opt_mem_per_job` | "60GB" | `null` |
| `intra_reaction_parallel.s1.crest_cores` | 20 | `null` |
| `intra_reaction_parallel.s1.sp_cores_per_job` | 20 | `null` |
| `intra_reaction_parallel.s1.sp_maxcore_per_job` | 5200 | `null` |
| `intra_reaction_parallel.s2.nproc` | 20 | `null` |
| `intra_reaction_parallel.s3.ts_opt_cores` | 20 | `null` |
| `intra_reaction_parallel.s3.ts_opt_mem` | "60GB" | `null` |
| `intra_reaction_parallel.s3.intermediate_opt_cores` | 20 | `null` |
| `intra_reaction_parallel.s3.intermediate_opt_mem` | "60GB" | `null` |
| `intra_reaction_parallel.s3.sp_cores_per_job` | 20 | `null` |
| `intra_reaction_parallel.s3.sp_maxcore_per_job` | 2600 | `null` |

#### 2.3 同步全部 6 个 YAML 文件

| 文件 | 修改方式 |
|------|---------|
| `defaults.yaml` | 全部 23+ 处 → null（主模板） |
| `defaults_alcl3_test.yaml` | 继承主模板，仅修改测试特有字段 |
| `defaults_la_test.yaml` | 继承主模板 |
| `defaults_la_test_dataset.yaml` | 继承主模板 |
| `defaults_la_test_fast.yaml` | 继承主模板 |
| `validation.yaml` | 继承主模板 |

**目标**：所有配置文件使用同一套 `null` 继承语义，不再各自维护 nproc 值。

#### 2.4 资源配置继承 + 预算分配代码实现

在 `rph_core/utils/resource_utils.py` 中新增 `resolve_resources(config)` 函数（**不新建模块**，避免两套解析逻辑）：

```python
def resolve_resources(config: dict) -> dict:
    """
    统一资源解析函数。在 orchestrator 启动时调用一次。
    
    执行：
    1. 若 resources.* 为 null → 自动检测系统 CPU/内存
    2. 各子区段（theory.*, step*.*, intra_reaction_parallel.*）若为 null → 继承 resources.*
    3. 并行 lane 做预算分配：per_job_cores = total_cores / max(workers, 1)
    4. 所有显式数值保持不动（override 语义）
    5. 写入 provenance
    """
```

**并行 lane 预算分配逻辑**（关键修正）：

```python
# intra_reaction_parallel 各 lane 的 *cores_per_job 继承自:
per_job_cores = total_cores / lane_workers

# 例如:
# total_cores=32, sp_workers=2 → sp_cores_per_job = 16
# total_cores=16, sp_workers=1 → sp_cores_per_job = 16
# total_cores=32, molecule_workers=2 → opt_cores_per_job = 16
```

**不把所有 lane 都填上全局核数**。`intra_reaction_scheduler.py` 中 `ResourceLane` 的默认值应改为 `None`，由 `resolve_resources()` 统一分配。

#### 2.5 集成到启动流程

```python
# orchestrator.py — ReactionProfileHunter.__init__()
self.config = load_config(config_path)

# Step 1: 统一资源解析（继承 + 预算 + 自动检测）
from rph_core.utils.resource_utils import resolve_resources
resolve_resources(self.config)

# Step 2: 可执行文件发现（已存在的 resolve_executable_config 扩展）
for key in _ALL_EXECUTABLES:
    resolve_executable_config(self.config, key, ..., allow_discovery=True)

# Step 3: 规范化（已有）
self.config, qc_fixes = normalize_qc_config(self.config, auto_fix=True)
```

---

### 阶段 3: 统一接口默认值

#### 3.1 原则

所有接口构造器参数默认值使用 `None`，在运行时由 `resolve_resources()` 解析填充。**不再添加新的硬编码默认值**。

| 接口 | 参数 | 旧默认值 | 新默认值 |
|------|------|---------|---------|
| `XTBInterface.__init__` | `nproc` | 1 | `None` |
| `CRESTInterface.__init__` | `nproc` | 1 | `None` |
| `XTBRunner.*` | nproc 回退 | `resources.get('nproc', 1)` | `resources.get('nproc', os.cpu_count())` |
| `GaussianInterface.__init__` | `nprocshared` | 16 | `None` |
| `GaussianInterface.__init__` | `mem` | "32GB" | `None` |
| `BernyTSDriver.__init__` | `nprocshared` | 16 | `None` |
| `QST2RescueDriver.__init__` | `nprocshared` | 16 | `None` |
| `IRCDriver.__init__` | `nprocshared` | 16 | `None` |
| `ResourceLane` | `nproc` | 8 | `None` |
| `ResourcePlan` | `qc_nproc` | 8 | `None` |
| `ResourcePlan` | `orca_maxcore` | 3000 | `None` |
| `ts_optimizer.py` | XTB nproc 回退 | 8 | `resources.nproc` 或 16 |

#### 3.2 调用点修改

每个构造器被调用时，调用方应传入解析后的值：

```python
# 旧
self.crest_interface = CRESTInterface(
    gfn_level=2, nproc=self.nproc, config=config)

# 新
self.crest_interface = CRESTInterface(
    gfn_level=2, nproc=config['step1']['crest']['threads'], config=config)
```

---

### 阶段 4: 扩展可执行文件发现层

#### 4.1 原则

**不新建模块**，在 `resource_utils.py` 中扩展 `resolve_executable_config()`，保持唯一解析入口。

#### 4.2 扩展 `resolve_executable_config()`

```python
def resolve_executable_config(
    config, program_key, env_vars=None,
    known_dirs=None, binary_names=None, fallback_paths=None,
    allow_discovery=True,
):
    """
    扩展后的解析顺序：
    1. config['executables'][key]['path']
       → 若文件存在，直接返回（显式路径优先）
       → 若文件不存在且 allow_discovery=True → 继续搜索
       → 若文件不存在且 allow_discovery=False → 输出 ERROR，fail-fast
    2. 检查环境变量（已有）
    3. shutil.which()（已有）
    4. 搜索 known_dirs（新增）
    5. 检查 fallback_paths（已有）
    6. 若全部失败 → WARNING, 返回 None
    """
```

**已知安装目录表**：

| program_key | known_dirs |
|------------|-----------|
| `gaussian` | `/opt/software/gaussian/*/g16`, `/opt/g16/g16`, `$HOME/g16/g16` |
| `orca` | `/opt/software/orca*/orca`, `/opt/orca/orca`, `$HOME/orca/orca` |
| `xtb` | `/opt/software/xtb*/bin/xtb`, `/opt/xtb/bin/xtb`, `$HOME/xtb/bin/xtb` |
| `crest` | `/opt/software/crest*/crest`, `/opt/crest/crest`, `$HOME/crest/crest` |
| `isostat` | `/opt/software/molclus*/isostat`, `/opt/molclus/isostat`, `$HOME/molclus*/isostat` |
| `shermo` | `/opt/software/shermo*/Shermo`, `/opt/shermo/Shermo`, `$HOME/shermo*/Shermo` |

#### 4.3 失效显式路径处理

如果 YAML 中的显式路径无效（文件不存在），不再静默覆盖：

```python
if exe_path and not Path(exe_path).is_file():
    if allow_discovery:
        logger.warning(
            f"Config path '{exe_path}' invalid. Auto-discovery enabled; "
            f"searching alternatives..."
        )
        # 写入 provenance: auto_discovered_from_invalid_path
    else:
        logger.error(
            f"Config path '{exe_path}' not found and "
            f"allow_discovery=False. Set the correct path manually."
        )
        # fail-fast
```

---

### 阶段 5: 验证

新增 `tests/test_config_unification.py`：

```python
def test_resolve_resources_inherits_nproc():
    """when theory.optimization.nproc=None, resources.nproc=32 → theory.optimization.nproc=32"""

def test_budget_allocation_with_workers():
    """total_cores=32, sp_workers=2 → sp_cores_per_job=16"""

def test_explicit_override_not_overwritten():
    """theory.optimization.nproc=8 → 不被 resources.nproc 覆盖"""

def test_invalid_path_discovery_does_not_silently_override():
    """executables.gaussian.path 指向不存在文件 → WARNING + 自动发现，日志记 provenance"""

def test_gau_xtb_writes_nproc():
    """GauXTBInterface(nproc=8).write_input_file() → %nprocshared=8"""

def test_all_config_nproc_null():
    """全部 6 个 YAML 文件的 theory.*.nproc 为 null"""
```

---

## 三、修改总表

### 3.1 Python 文件修改

| 文件 | 修改 | 阶段 |
|------|------|------|
| `gau_xtb_interface.py` | `%nprocshared=1` → `self.nproc` | 1.1 |
| `runners.py` | `nproc=1` → `resources.nproc` | 1.1 |
| `resource_utils.py` | `safety_factor=0.8→0.65`, 新增 `resolve_resources()`, 扩展 `resolve_executable_config()` | 1.2, 2.4, 4.2 |
| `orca_interface.py` | `0.8→0.65` | 1.2 |
| `engine.py` | `'16GB'→'32GB'`（4 处） | 1.3 |
| `qc_interface.py` | `'48GB'→'32GB'`, XTBInterface/CRESTInterface nproc 默认 `None` | 1.4, 3 |
| `orchestrator.py` | 集成 `resolve_resources()` + `resolve_executable_config()` 扩展 | 2.5 |
| `xtb_runner.py` | `1→os.cpu_count()` | 3 |
| `ts_optimizer.py` | XTB nproc 回退统一 | 3 |
| `intra_reaction_scheduler.py` | ResourceLane 默认 `None` | 3 |
| `scheduler.py` | ResourcePlan 默认 `None` | 3 |
| `berny_driver.py` | nprocshared/mem 默认 `None` | 3 |
| `qst2_rescue.py` | nprocshared/mem 默认 `None` | 3 |
| `irc_driver.py` | nprocshared/mem 默认 `None` | 3 |

### 3.2 YAML 配置文件修改

| 文件 | 修改 | 阶段 |
|------|------|------|
| `defaults.yaml` | 23+ 处 nproc/mem/threads → `null` | 2.1 |
| `defaults_alcl3_test.yaml` | 同步 | 2.2 |
| `defaults_la_test.yaml` | 同步 | 2.2 |
| `defaults_la_test_dataset.yaml` | 同步 | 2.2 |
| `defaults_la_test_fast.yaml` | 同步 | 2.2 |
| `validation.yaml` | 同步 | 2.2 |

---

## 四、设计决策记录

| 决策 | 选项 | 选择 | 理由 |
|------|------|------|------|
| 自动发现放哪里 | 新建模块 vs 扩展现有 | **扩展 `resource_utils.resolve_executable_config()`** | 避免两套解析逻辑、两套日志 |
| 接口默认值 | 1/8/16 vs None | **`None`** | 不再制造新的硬编码 |
| 失效路径处理 | 静默覆盖 vs fail-fast | **WARNING + 可选 discovery + provenance** | 平衡可用性与安全性 |
| 配置文件修改 | 逐个 vs 全部 | **全部 6 个一起** | 避免配置分叉 |
| 修改顺序 | 大重构 vs 增量 | **先修 BUG（阶段1），再统一（阶段2-4）** | 每步可测，风险可控 |
| 并行 lane 分配 | 全继承 vs 预算分配 | **预算分配（total_cores / workers）** | 防止超订阅 |

