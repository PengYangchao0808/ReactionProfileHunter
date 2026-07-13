# RPH 配置修复 — 后续阶段计划

> **版本**: v3.0 | **日期**: 2026-06-30 | **状态**: 待实施

---

## 总体路线

```
当前状况 → Phase 1 [P1] 资源闭环 → Phase 2 [P3] 软件发现 → Phase 3 [P2] LA切换 → Phase 4 [P2] 测试
  (75%完成)     (3 个残留, ~15行)    (~90行新代码)          (~50行)                (~120行)
```

**顺序执行原则**：Phase 1 未完成前不进入 Phase 2/3（否则 S3 并行仍可能超配）。

---

## Phase 1：资源调度闭环修复 [P1]

**目标**：使 `resolve_resources()` 成为真正的资源单一入口。

### 1.1 返回值修正

**问题**：`resource_utils.py:303` 的 `resolve_resources()` 函数签名声明 `-> dict`，但函数末尾没有 `return config`，实际返回 `None`。后续测试/脚本复用会踩坑。

**修改**：

```python
# resource_utils.py — resolve_resources() 末尾
    _resolve_irp_defaults(irp, resources, nproc, mem)
    return config   # ← 新增，与签名一致
```

**涉及文件**：`rph_core/utils/resource_utils.py:407`

---

### 1.2 `_resolve_irp_defaults()` 回写 total_cores/total_mem

**问题**：`intra_reaction_scheduler.py:113` 读取 `irp.get('total_cores', 16)` 和 `irp.get('total_mem', '52GB')`，但 `resolve_resources()` 从未回写这两个字段，导致预算检查使用旧默认值 16 cores / 52GB，与 `resources.mem=32GB` 不一致。

**修改**：

```python
# resource_utils.py — _resolve_irp_defaults() 末尾新增
    irp['total_cores'] = total_cores
    irp['total_mem'] = total_mem_str
```

**涉及文件**：`rph_core/utils/resource_utils.py:355`（函数末尾）

---

### 1.3 S3 TS/Intermediate 并行语义修正

**问题**：`ts_optimizer.py:742` 在 `parallel_intermediate_ts=true` 时会并行提交 TS 和 Intermediate。当前 `_resolve_irp_defaults()` 给两者都分配 full nproc，导致超配。

**规则**：

```
parallel_intermediate_ts = false（串行）:
  ts_cores = resources.nproc
  inter_cores = resources.nproc
  → 从不同时运行，安全

parallel_intermediate_ts = true（并行）:
  ts_cores : inter_cores = 3 : 1（历史经验值）
  例: 16核 → ts=12, inter=4
  内存按相同比例拆分
```

**修改 1** — `resource_utils.py:_resolve_irp_defaults()`：

```python
# S3 ──
parallel_it = s3.get('parallel_intermediate_ts', True)

if s3.get('ts_opt_cores') is None:
    if parallel_it:
        s3['ts_opt_cores'] = max(1, int(total_cores * 3 / 4))   # 3/4 给 TS
    else:
        s3['ts_opt_cores'] = total_cores

if s3.get('intermediate_opt_cores') is None:
    if parallel_it:
        s3['intermediate_opt_cores'] = max(1, int(total_cores / 4))  # 1/4 给 Inter
    else:
        s3['intermediate_opt_cores'] = total_cores

# 内存必须按相同比例拆分，避免并行时总内存超 budget
if s3.get('ts_opt_mem') is None:
    if parallel_it:
        ts_mem_gb = max(8, int(mem_to_mb(total_mem_str) * 3 / 4 / 1024))
        s3['ts_opt_mem'] = f"{ts_mem_gb}GB"
    else:
        s3['ts_opt_mem'] = total_mem_str
if s3.get('intermediate_opt_mem') is None:
    if parallel_it:
        inter_mem_gb = max(4, int(mem_to_mb(total_mem_str) / 4 / 1024))
        s3['intermediate_opt_mem'] = f"{inter_mem_gb}GB"
    else:
        s3['intermediate_opt_mem'] = total_mem_str
```

**修改 2** — `ts_optimizer.py:511`：

```python
# 将硬编码日志 "TS (12 cores) + Intermediate (4 cores)" 改为动态值
self.logger.info(
    f"TS ({ts_lane.nproc} cores) + Intermediate ({int_lane.nproc} cores)"
)
```

**涉及文件**：`rph_core/utils/resource_utils.py`，`rph_core/steps/step3_opt/ts_optimizer.py`

---

### 1.4 Scheduler overcommit 行为可配置

**问题**：`intra_reaction_scheduler.py:231` 超配时仅 warning 后继续执行。

**方案**：新增配置开关和对应逻辑

```yaml
intra_reaction_parallel:
  on_overcommit: "warn"   # warn | fail | fallback_sequential
```

代码逻辑：

```python
# intra_reaction_scheduler.py:231
overcommit_policy = parallel_cfg.get('on_overcommit', 'warn')
if total_requested > total_cores or total_mem_requested > total_mem_mb:
    if overcommit_policy == 'fail':
        raise RuntimeError(f"Resource overcommit detected: ...")
    elif overcommit_policy == 'fallback_sequential':
        logger.warning("Overcommit: falling back to sequential execution")
        # 内部直接降级为串行，不改变调用方接口
        return self._run_sequential(jobs)
    else:  # 'warn'（默认）
        logger.warning(f"Resource overcommit: {total_requested} > {total_cores} cores")
        # 继续当前并行路径
```

**涉及文件**：`rph_core/utils/intra_reaction_scheduler.py`

---

### 1.5 最终预算校验 `_validate_irp_budget()`

**问题**：当前方案主要修 `null` 自动推导。但如果用户**显式**写 `ts_opt_cores: 16`、`intermediate_opt_cores: 16`，仍会超配。warn/fail/fallback_sequential 应基于合并 override 后的最终状态，而非仅自动生成值。

**方案**：在 `_resolve_irp_defaults()` 末尾统一校验：

```python
def _resolve_irp_defaults(irp, resources, nproc, mem):
    # ... 前面的默认值填充逻辑 ...

    # ── 最终预算校验（覆盖自动推导和显式 override 的最终值）
    _validate_irp_budget(irp, total_cores, total_mem_mb)
```

校验函数：

```python
def _validate_irp_budget(irp, total_cores, total_mem_mb):
    """校验所有 lane 合并后不超总预算。"""
    s3 = irp.get('s3', {})
    parallel_it = s3.get('parallel_intermediate_ts', True)

    # S3 并行时 TS+Inter 不得超 total_cores
    if parallel_it:
        ts_cores = s3.get('ts_opt_cores', 0) or 0
        inter_cores = s3.get('intermediate_opt_cores', 0) or 0
        if ts_cores + inter_cores > total_cores:
            logger.warning(
                f"S3 budget exceeded: TS({ts_cores})+Inter({inter_cores})="
                f"{ts_cores+inter_cores} > total({total_cores})"
            )

    # SP 总内存校验
    for lane_name in ('s1', 's3'):
        lane = irp.get(lane_name, {})
        sp_cores = lane.get('sp_cores_per_job', 0) or 0
        sp_maxcore = lane.get('sp_maxcore_per_job', 0) or 0
        sp_workers = lane.get('sp_workers', 1)
        total_sp_mem = sp_cores * sp_maxcore * sp_workers
        budget_mem = int(total_mem_mb * irp.get('resources', {}).get('orca_maxcore_safety', 0.65))
        if total_sp_mem > budget_mem:
            logger.warning(
                f"{lane_name} SP memory {total_sp_mem/1024:.1f}GB > "
                f"budget {budget_mem/1024:.1f}GB"
            )
```

**涉及文件**：`rph_core/utils/resource_utils.py`（新增函数，~25 行）

---

### 1.6 GauXTBOptimizer 继承 resources.nproc

**问题**：`gau_xtb_interface.py:414` 的 `GauXTBOptimizer.__init__` 默认 `nproc=1`，外部直接调用 `GauXTBOptimizer(config)` 时回到单核。

**修改**：

```python
# gau_xtb_interface.py — GauXTBOptimizer.__init__
def __init__(self, config, nproc=None, ...):
    if nproc is None:
        nproc = int(config.get('resources', {}).get('nproc', os.cpu_count() or 1))
    self.nproc = nproc
```

与已修复的 `GauXTBInterface.__init__`（line 67）保持相同模式。

**涉及文件**：`rph_core/utils/gau_xtb_interface.py`

---

## Phase 2：软件自动发现 [P3]

**目标**：在现有 `resource_utils.resolve_executable_config()` 上扩展，无需新建平行体系。

### 2.1 新增配置段

```yaml
# config/defaults.yaml 新增
executables:
  gaussian: { path: "/opt/software/gaussian/g16/g16" }
  orca:     { path: "/opt/software/orca/orca" }
  xtb:      { path: "/opt/software/xtb/bin/xtb" }
  crest:    { path: "/opt/software/crest/crest" }
  isostat:  { path: "/opt/software/molclus/isostat" }
  shermo:   { path: "/opt/software/shermo/Shermo" }

  discovery:
    enabled: true                      # 是否启用自动发现
    fail_on_invalid_explicit: true     # 显式路径无效时是否报错
    known_dirs:
      gaussian: ["/opt/software/gaussian", "/opt/g16", "/usr/local/g16"]
      orca:     ["/opt/software/orca", "/opt/orca"]
      xtb:      ["/opt/software/xtb", "/opt/xtb", "/usr/local/bin"]
      crest:    ["/opt/software/crest", "/opt/crest", "/usr/local/bin"]
      isostat:  ["/opt/software/molclus", "/opt/molclus"]
      shermo:   ["/opt/software/shermo", "/opt/shermo"]
      mpirun:   ["/opt/openmpi", "/usr/lib/openmpi"]
```

### 2.2 扩展 `resolve_executable_config()`

**位置**：`resource_utils.py:189`

**解析优先级**（新增第 4 级，source 返回值明确化）：

```
1. YAML config['executables'][key]['path']
   → 若文件存在且可执行 → 返回，source="explicit"
   → 若文件不存在且 fail_on_invalid_explicit=true → 抛出 RuntimeError
   → 若文件不存在且 fail_on_invalid_explicit=false → 继续搜索（降级模式）
2. 环境变量（已有）→ source="env"
3. shutil.which()（已有）→ source="path"
4. known_dirs glob 搜索（新增）→ source="known_dirs"
   → 遍历预设目录列表，用 Path.glob 匹配 binary_name
5. 全部失败 → 明确报错，列出搜索路径，source="not_found"
```

### 2.3 新增 `resolve_all_executables()`

```python
# resource_utils.py 新增
def resolve_all_executables(config) -> Dict[str, ResolveResult]:
    """一次性解析全部 6+1 个可执行文件。

    在 orchestrator 启动时调用一次，位于 resolve_resources() 之后。

    Returns:
        Dict[program_key, ResolveResult]
        ResolveResult = {path: Path|None, source: str, found: bool}
    """
    programs = ['gaussian', 'orca', 'xtb', 'crest', 'isostat', 'shermo', 'mpirun']
    results = {}
    for prog in programs:
        results[prog] = resolve_executable_config(
            config, prog,
            env_vars=_ENV_MAP[prog],
            known_dirs=_KNOWN_DIRS.get(prog, []),
            binary_names=_BINARY_NAMES[prog],
            fail_on_invalid=config.get('executables', {}).get('discovery', {}).get('fail_on_invalid_explicit', True),
        )
    return results
```

### 2.4 集成到 orchestrator

```python
# orchestrator.py:145-153
self.config = load_config(config_path)
resolve_resources(self.config)              # Phase 1: 资源（已有）

# Phase 2: 可执行文件发现（新增）
from rph_core.utils.resource_utils import resolve_all_executables
exe_results = resolve_all_executables(self.config)
for prog, result in exe_results.items():
    if result.found:
        logger.info(f"[discovery] {prog}: {result.path} ({result.source})")
    else:
        logger.warning(f"[discovery] {prog}: NOT FOUND")

self.config, qc_fixes = normalize_qc_config(self.config, auto_fix=True)
```

### 2.5 Resolution 报告输出示例

```
[discovery] Gaussian:   /opt/software/gaussian/g16/g16  (explicit)
[discovery] ORCA:      /opt/software/orca/orca          (explicit)
[discovery] xTB:       /opt/software/xtb/bin/xtb        (path)
[discovery] CREST:     /opt/software/crest/crest         (known_dirs)
[discovery] ISOSTAT:   /opt/software/molclus/isostat     (known_dirs)
[discovery] Shermo:    /opt/software/shermo/Shermo       (known_dirs)
[discovery] MPI:       /opt/openmpi418/bin/mpirun        (path)
```

### 2.6 required vs optional 可执行文件区分

**问题**：gaussian/orca/xtb/crest/isostat/shermo/mpirun 不是所有 run 都必需。若 `resolve_all_executables()` 对所有缺失项 fail-fast，会阻断仅跑 xTB 或跳过某些阶段的任务。

**方案**：根据当前 config 的运行时上下文判断必需性：

```python
def _is_executable_required(config, program_key) -> bool:
    """根据配置判断该可执行文件是否必需。"""
    theory = config.get('theory', {})
    step1 = config.get('step1', {})

    REQUIRED_MAP = {
        'gaussian': lambda: theory.get('optimization', {}).get('engine') == 'gaussian'
                         or theory.get('single_point', {}).get('engine') == 'gaussian',
        'orca':     lambda: theory.get('optimization', {}).get('engine') == 'orca'
                         or theory.get('single_point', {}).get('engine') == 'orca',
        'xtb':      lambda: True,  # CREST/xTB scan 均依赖
        'crest':    lambda: not step1.get('conformer_search', {}).get('two_stage_enabled', False),
        'isostat':  lambda: True,  # CREST 聚类必需
        'shermo':   lambda: config.get('thermo', {}).get('engine', 'shermo') == 'shermo',
        'mpirun':   lambda: theory.get('single_point', {}).get('engine') == 'orca'
                         and (theory.get('single_point', {}).get('nproc') or 1) > 1,
    }
    return REQUIRED_MAP.get(program_key, lambda: False)()

# 在 resolve_all_executables() 中，仅必需项缺失时 error；非必需项缺失时仅 warning
if not result.found:
    if _is_executable_required(config, prog):
        logger.error(f"[discovery] REQUIRED {prog}: NOT FOUND")
    else:
        logger.warning(f"[discovery] optional {prog}: NOT FOUND (skipped)")
```

---

## Phase 3：强弱 LA 切换补全 [P2]

**目标**：从"能改 surrogate 名称"升级为"每个 surrogate 有独立参数"。

### 3.1 注册 BF₃

在 `config/defaults.yaml` 和所有测试配置的 `lewis_acid.surrogate_registry` 中添加：

```yaml
lewis_acid:
  surrogate_registry:
    LiCl:  { name: LiCl,  smiles: '[Li]Cl',       charge: 0, multiplicity: 1 }
    MgCl2: { name: MgCl2, smiles: 'Cl[Mg]Cl',      charge: 0, multiplicity: 1 }
    AlCl3: { name: AlCl3, smiles: 'Cl[Al](Cl)Cl',  charge: 0, multiplicity: 1 }
    BF3:   { name: BF3,   smiles: 'B(F)(F)F',      charge: 0, multiplicity: 1 }
```

**当前方案**：直接写入各 YAML 的 `surrogate_registry` 段。由于 `config_loader.py` 没有 `!include` 机制，先不抽取独立文件，避免新增加载器改造任务。

**更新测试断言**：

```python
# tests/test_lewis_acid_phase1.py:220
# 旧: set(registry.keys()) == {"LiCl", "MgCl2", "AlCl3"}
# 新: set(registry.keys()) == {"LiCl", "MgCl2", "AlCl3", "BF3"}
```

### 3.2 新增 `surrogate_overrides` 机制

```yaml
# config/defaults.yaml 新增
lewis_acid:
  enabled: true
  surrogate: AlCl3           # ← 切换点
  surrogate_registry: {...}

  # 基础质量门控（LiCl 默认）
  s1_sampling: { ... }
  s2_quality_gate: { ... }
  s3_quality_gate: { ... }

  # 各 surrogate 专属覆盖（新增）
  surrogate_overrides:
    BF3:
      s1_sampling:
        max_metal_substrate_distance: 3.5
      s2_quality_gate:
        metal_o_max: 2.3
        metal_cl_max: 3.0
      s3_quality_gate:
        additive_mode_projection_max: 0.30
    MgCl2:
      s2_quality_gate:
        metal_o_max: 3.0
        coordination_switch_allowed: true
      s3_quality_gate:
        additive_mode_projection_max: 0.20
```

### 3.3 `_prepare_lewis_acid_system()` deep merge

现有函数是 `orchestrator.py:196` 的**实例方法** `ReactionProfileHunter._prepare_lewis_acid_system(product_smiles, precursor_smiles)`，不是独立函数。deep merge 应在其内部处理 `self.config["lewis_acid"]`：

```python
# orchestrator.py — ReactionProfileHunter._prepare_lewis_acid_system() 内部
def _prepare_lewis_acid_system(self, product_smiles, precursor_smiles):
    from copy import deepcopy
    la_config = deepcopy(self.config.get('lewis_acid', {}))
    surrogate = la_config.get('surrogate', 'LiCl')

    # 1. 从 registry 加载基础定义（已有逻辑）
    registry = la_config.get('surrogate_registry', {})
    surrogate_def = registry.get(surrogate, {})

    # 2. 应用 surrogate_overrides（新增：深合并到 la_config）
    overrides = la_config.get('surrogate_overrides', {}).get(surrogate, {})
    _deep_merge(la_config, overrides)
    # 覆盖 s1_sampling / s2_quality_gate / s3_quality_gate 中对应的字段

    # 3. 后续构建 LewisAcidAdditive（已有逻辑，使用合并后的 la_config）
    ...

### 3.4 命名建议

`LewisAcidAdditive` 中新增 `center_element`/`center_atom_index`，与 `metal_element`/`metal_atom_index` 做兼容别名：

```python
@property
def center_element(self) -> Optional[str]:
    """中心元素（metal_element 的通用别名）。"""
    return self.metal_element
```

---

## Phase 4：测试与验收 [P2]

### 4.1 新增测试文件

| 文件 | 测试内容 |
|------|---------|
| `tests/test_resource_resolution.py`（新建） | `resolve_resources()` 返回 config；`total_cores`/`total_mem` 回写；S3 parallel lane 总核数/内存不超过预算；SP maxcore worker 预算正确 |
| `tests/test_executable_resolution.py`（新建） | 显式路径优先；无效显式路径 fail-fast；PATH fallback；`known_dirs` fake executable discovery |
| `tests/test_lewis_acid_phase1.py`（扩展） | registry 包含 LiCl/MgCl2/AlCl3/BF3；BF3 的 B 为 center、F 为 halide；forming bonds 不包含 LA 原子；`surrogate_overrides` 生效 |

### 4.2 验收命令

```bash
# Import 风格门
python scripts/ci/check_imports.py rph_core

# 资源测试
pytest tests/test_resource_resolution.py -v

# 可执行文件发现测试
pytest tests/test_executable_resolution.py -v

# LA 测试
# LA 测试（需要 RDKit 环境）
pytest tests/test_lewis_acid_phase1.py -v -p no:rdkit 2>/dev/null || \
  conda run -n rdkit_env pytest tests/test_lewis_acid_phase1.py -v

# 全量（跳过 rdkit 依赖）
pytest tests/ --ignore=tests/test_mrrho.py -v
```

---

## 修改总表

### 本次修改（Phase 1）

| 修改 | 文件 | 改动量 | 说明 |
|------|------|--------|------|
| `return config` | `resource_utils.py:resolve_resources` | +1 行 | 与签名一致 |
| 回写 total_cores/total_mem | `resource_utils.py:_resolve_irp_defaults` | +2 行 | 使 budget 检查使用正确值 |
| S3 parallel 3:1 拆分 | `resource_utils.py:_resolve_irp_defaults` | ~10 行 | TS=12, Inter=4 (16核时) |
| 日志改动态 | `ts_optimizer.py:511` | ~2 行 | 硬编码日志改动态 |
| overcommit 可配置 | `intra_reaction_scheduler.py:231` | ~15 行 | warn/fail/fallback_sequential |
| GauXTBOptimizer 继承 | `gau_xtb_interface.py:414` | ~5 行 | nproc=None 自动继承 |
| **Phase 1 合计** | **4 文件** | **~35 行** | |

### 本次修改（Phase 2-4，后续实施）

| Phase | 文件数 | 代码新增 | 主要风险 |
|-------|--------|---------|---------|
| 2 软件发现 | 3 | ~90 行 | 中 — 新功能，需测试 |
| 3 LA 切换 | 8 | ~50 行 | 中 — 配置分散 |
| 4 测试 | 3 | ~120 行 | 低 — 纯新增 |

---

## 当前最小建议

**先做 Phase 1**（~1 小时工作量）。完成后确认：

```bash
python -c "
from rph_core.utils.resource_utils import resolve_resources
import yaml
with open('config/defaults.yaml') as f:
    cfg = yaml.safe_load(f)
result = resolve_resources(cfg)
assert result is not None, 'resolve_resources() must return config'
print('✅ Phase 1 闭环验证通过')
"
```

Phase 1 闭环后，"资源自动调度"才算真正完成。之后再进入 Phase 2（软件发现）和 Phase 3（LA 切换）。
