# S2 内嵌式 fast-3c TS 预优化重构方案

## 1. 纠偏说明

本方案修正一个关键理解错误：

- **复用 S3** 不是指 S2 预优化必须使用 S3 最终 DFT 理论层，例如 `M062X/def2-SVP`。
- **复用 S3** 指复用 S3 已有的优化调度、结果解析、XYZ 写出、失败分型、fallback、metadata 和测试框架。
- S2 内嵌预优化的理论层仍应是默认快速 3c 方法，例如 `r2SCAN-3c` 或 `B97-3c`。

因此，正确目标是：

```text
S2 内部 preopt
  = fast-3c 方法
  + S3-style optimization runner / result handling
  + TS guess 使用 TS-like 优化关键词
  + intermediate 使用 normal opt
  + 不新增 S2.5 数据流
```

当前测试暴露的问题不是“3c 方法不该用”，而是当前代码把 TS guess 送进了：

```text
! r2SCAN-3c Opt ...
%geom Constraints ...
```

这只是普通约束基态优化，不是 TS 预优化。它和 TS 搜索不沾边，所以测试结果没有意义。

## 2. 总体目标

1. 不新增 `S2.5` 阶段。
2. 不新增 `run_step2_preopt()`。
3. 不改 orchestrator 的 S2/S3 阶段边界。
4. 预优化完全内嵌在 `run_step2()` 内，属于 S2 refinement 子步骤。
5. S2 对外仍只返回 `Step2Artifacts`。
6. 预优化理论层默认使用 fast-3c：
   - 首选 `r2SCAN-3c`
   - 可选 `B97-3c`
7. 最大化复用 S3 的计算执行框架，而不是手写一套孤立 ORCA 调用。
8. TS guess 预优化必须使用 TS-like 关键词，不允许再用普通 `Opt` 冒充。

## 3. 当前实现的问题

当前 `S2PreOptDriver` 的核心问题：

```python
orca = self._create_orca_interface()
result = orca.constrained_optimize(
    xyz_file=input_xyz,
    output_dir=output_dir,
    charge=charge,
    spin=spin,
    constraints_block=constraints_block,
    timeout=self.timeout,
)
```

而 `ORCAInterface.constrained_optimize()` 内部固定使用：

```python
route = self._render_route(task_type="opt")
```

生成的是：

```text
! r2SCAN-3c Opt tightSCF TightOpt ... CPCM(Acetone)
```

这会导致：

- TS guess 被当作普通最小值优化处理。
- forming bonds 虽然被冻结，但 Hessian/TS 搜索逻辑完全没有进入。
- `preopt_meta.json` 记录的是 fast-3c normal constrained opt，而不是 TS preopt。
- S3 后续失败/成功无法归因于“TS 初猜改善”。

## 4. 正确设计

### 4.1 数据流仍然属于 S2

```text
run_step2()
  ├─ xTB retro_scan
  ├─ optional xTB path_search
  ├─ raw S2 outputs
  │    ├─ ts_guess.xyz
  │    └─ intermediate.xyz
  ├─ embedded fast-3c preopt
  │    ├─ backup raw S2 outputs
  │    ├─ TS guess: fast-3c TS-like preopt
  │    ├─ intermediate: fast-3c normal opt
  │    ├─ accept or fallback raw
  │    └─ write preopt_meta.json
  ├─ write S2 provenance
  └─ return Step2Artifacts
```

对外不新增阶段，S3 仍然只消费：

```text
S2_Retro/ts_guess.xyz
S2_Retro/intermediate.xyz
```

### 4.2 复用 S3 的含义

应复用这些 S3 基础设施：

- `QCTaskRunner` 的统一优化结果模型。
- `QCOptimizationResult`。
- `_write_optimized_xyz()` 原子序与符号链保护。
- `_count_imaginary()` 虚频统计。
- `_coalesce_qc_error()` 错误信息收集。
- `_analyze_log_for_fatal_errors()` 致命错误判断。
- S3 rescue/raw fallback 的思路和 metadata 风格。
- ORCA/Gaussian interface 的已有执行、解析、sandbox、timeout 能力。

不应复用 S3 的最终 DFT method/basis。

### 4.3 fast-3c backend

新增一个专门的 fast-3c preopt 配置覆盖层，构造临时 runner/config：

```yaml
method: r2SCAN-3c
basis: ""
aux_basis: ""
engine: orca
solvent: acetone
solvent_model: CPCM
```

也可以切换：

```yaml
method: B97-3c
```

该配置只服务 S2 内嵌 preopt，不污染 S3 最终 DFT 配置。

## 5. 推荐配置

替换现有 `step2.pre_opt` 为以下语义：

```yaml
step2:
  pre_opt:
    enabled: false
    backend: fast_3c_s3_style
    engine: orca
    method: r2SCAN-3c          # r2SCAN-3c | B97-3c
    basis: ""
    aux_basis: ""
    solvent: acetone
    solvent_model: CPCM
    nproc: null
    charge: null
    multiplicity: null
    timeout_seconds: 7200
    fallback_to_raw_on_failure: true

    trigger:
      always: false
      on_lewis_acid: true
      on_status: [DEGRADED, FAILED]

    ts_guess:
      enabled: true
      mode: optts
      route_task: ts_freq       # ts | ts_freq
      max_cycles: 50
      calc_hess: true
      recalc_hess_every: null
      freeze_forming_bonds: false
      require_one_imaginary: false
      accept_zero_imaginary: false
      fallback_to_raw_on_failure: true

    intermediate:
      enabled: true
      mode: normal_opt
      route_task: opt           # opt | opt_freq
      max_cycles: 100
      fallback_to_raw_on_failure: true
```

说明：

- `ts_guess.mode: optts` 是默认推荐。TS guess 必须进入 TS 优化语义。
- `freeze_forming_bonds` 默认建议先设为 `false`。因为 TS 优化本身需要沿反应坐标找负曲率，冻结 forming bonds 可能阻断 TS 模式。
- 若后续证明确实需要约束，可再实现 ORCA `OptTS + Constraints` 的组合，但不能回退到普通 `Opt`。
- `intermediate` 是基态结构，可使用普通 `Opt`。

## 6. ORCA 关键词要求

### 6.1 TS guess preopt

目标输入应类似：

```text
! r2SCAN-3c OptTS Freq tightSCF TightOpt noautostart miniprint nopop CPCM(Acetone)
%geom
  MaxIter 50
  Calc_Hess true
end
* xyz 0 1
...
*
```

或 B97-3c：

```text
! B97-3c OptTS Freq tightSCF TightOpt noautostart miniprint nopop CPCM(Acetone)
```

禁止再生成：

```text
! r2SCAN-3c Opt ...
```

用于 `ts_guess`。

### 6.2 intermediate preopt

目标输入可为：

```text
! r2SCAN-3c Opt tightSCF TightOpt noautostart miniprint nopop CPCM(Acetone)
%geom
  MaxIter 100
end
* xyz 0 1
...
*
```

intermediate 是 normal opt，此处使用 `Opt` 是合理的。

## 7. 代码修改方案

### 7.1 不再直接在 S2PreOptDriver 中裸调 ORCA

文件：

```text
rph_core/steps/step2_retro/preopt_driver.py
```

删除或弃用：

```python
_create_orca_interface()
_build_geom_block()
orca.constrained_optimize(...)
```

替换为：

```python
from rph_core.utils.qc_task_runner import QCTaskRunner

runner = QCTaskRunner.from_config(self._build_fast3c_runner_config())
```

然后调用新增的 preopt runner API：

```python
ts_result = runner.run_fast3c_ts_preopt_cycle(...)
int_result = runner.run_fast3c_normal_preopt_cycle(...)
```

### 7.2 在 QCTaskRunner 中新增 fast-3c preopt API

文件：

```text
rph_core/utils/qc_task_runner.py
```

新增：

```python
def run_fast3c_ts_preopt_cycle(
    self,
    xyz_file: Path,
    output_dir: Path,
    charge: int = 0,
    spin: int = 1,
    max_cycles: int = 50,
    include_freq: bool = True,
    timeout: Optional[int] = None,
) -> QCOptimizationResult:
    ...
```

职责：

1. 要求 `engine_type == "orca"`。
2. 通过 ORCAInterface 生成 `OptTS` route。
3. 注入 `%geom MaxIter`。
4. 可选注入 `Calc_Hess true`。
5. 执行 ORCA TS preopt。
6. 解析 energy、coordinates、frequencies。
7. 使用 `_write_optimized_xyz()` 写出 XYZ。
8. 统计 imaginary_count。
9. 返回 `QCOptimizationResult(method_used="Fast3C_OptTS_PreOpt")`。

### 7.3 ORCAInterface 需要新增通用 optimize_with_task

文件：

```text
rph_core/utils/orca_interface.py
```

当前问题是：

```python
constrained_optimize() -> task_type="opt"
_run_normal_optimization() -> task_type="opt_freq"
ts_optimization() -> 独立路径，但不方便传 MaxIter/fast preopt 策略
```

建议新增一个通用方法：

```python
def optimize_with_task(
    self,
    xyz_file: Path,
    output_dir: Path,
    task_type: str,
    charge: int,
    spin: int,
    geom_block: Optional[str] = None,
    timeout: Optional[int] = None,
) -> QCResult:
    ...
```

其中 `task_type` 允许：

```text
opt
opt_freq
ts
ts_freq
```

内部直接复用：

```python
route = self._render_route(task_type=task_type)
cpcm_block = self._render_cpcm_block()
rendered_blocks = self._render_named_blocks(self.render_blocks())
```

这样 TS preopt 可调用：

```python
orca.optimize_with_task(
    xyz_file=ts_guess_raw,
    output_dir=preopt_dir / "ts_guess",
    task_type="ts_freq",
    charge=charge,
    spin=spin,
    geom_block="%geom\n  MaxIter 50\n  Calc_Hess true\nend",
    timeout=timeout,
)
```

intermediate preopt 可调用：

```python
orca.optimize_with_task(
    xyz_file=intermediate_raw,
    output_dir=preopt_dir / "intermediate",
    task_type="opt",
    charge=charge,
    spin=spin,
    geom_block="%geom\n  MaxIter 100\nend",
    timeout=timeout,
)
```

### 7.4 不推荐继续使用 constrained_optimize 处理 TS guess

`constrained_optimize()` 当前固定：

```python
task_type="opt"
```

因此它只能用于 normal constrained opt。除非改造成支持 `task_type="ts"`，否则不应用于 `ts_guess`。

如果确实要约束 forming bonds，应该新增：

```python
def constrained_optimize_with_task(..., task_type: str)
```

并明确生成：

```text
! r2SCAN-3c OptTS ...
%geom
  Constraints
    { B i j C }
  end
end
```

而不是：

```text
! r2SCAN-3c Opt ...
```

## 8. S2PreOptDriver 新职责

`S2PreOptDriver` 只负责 S2 语义，不负责底层 QC 细节：

1. 判断 trigger。
2. 备份 raw 文件。
3. 构造 fast-3c runner config。
4. 调用 QCTaskRunner preopt API。
5. 接受/拒绝优化结果。
6. 写 `preopt_meta.json`。
7. 返回更新后的 `ts_guess_xyz` 和 `intermediate_xyz`。

示意：

```python
def _build_fast3c_runner_config(self) -> Dict[str, Any]:
    cfg = copy.deepcopy(self.config)
    pre = self.pre_opt_cfg
    cfg["theory"]["optimization"] = {
        "engine": pre.get("engine", "orca"),
        "method": pre.get("method", "r2SCAN-3c"),
        "basis": pre.get("basis", ""),
        "aux_basis": pre.get("aux_basis", ""),
        "solvent": pre.get("solvent", "acetone"),
        "solvent_model": pre.get("solvent_model", "CPCM"),
    }
    return cfg
```

注意：

- 这是临时 config，只用于 S2 preopt。
- 不应修改 `hunter.config` 原对象。
- 不应影响 S3 的 `theory.optimization`。

## 9. 接受/回退策略

### 9.1 TS guess

默认推荐：

- ORCA 正常结束且坐标完整：可接受。
- 如果 `include_freq=true`：
  - `imaginary_count == 1`：接受。
  - `imaginary_count == 0`：默认拒绝并回退 raw。
  - `imaginary_count > 1`：默认可配置，建议先拒绝。
- ORCA error、无坐标、原子数不一致：回退 raw。

配置：

```yaml
require_one_imaginary: false
accept_zero_imaginary: false
```

解释：

- `require_one_imaginary=false` 是因为 fast-3c preopt 不是最终 TS 验证。
- `accept_zero_imaginary=false` 是为了防止 preopt 已经坍塌到 minimum 还覆盖 raw TS guess。

### 9.2 intermediate

- normal opt 收敛且坐标完整：接受。
- 失败则回退 raw。

## 10. preopt_meta.json schema

路径：

```text
S2_Retro/preopt/preopt_meta.json
```

建议：

```json
{
  "schema_version": "s2_embedded_fast3c_preopt_v1",
  "embedded_in_step": "s2",
  "backend": "fast_3c_s3_style",
  "engine": "orca",
  "method": "r2SCAN-3c",
  "basis": "",
  "aux_basis": "",
  "solvent": "acetone",
  "solvent_model": "CPCM",
  "trigger_reason": "lewis_acid_present",
  "charge": 0,
  "multiplicity": 1,
  "raw": {
    "ts_guess_xyz": "S2_Retro/ts_guess_raw_s2.xyz",
    "intermediate_xyz": "S2_Retro/intermediate_raw_s2.xyz"
  },
  "ts_guess": {
    "attempted": true,
    "accepted": true,
    "method_used": "Fast3C_OptTS_PreOpt",
    "task_type": "ts_freq",
    "route": "! r2SCAN-3c OptTS Freq tightSCF TightOpt CPCM(Acetone)",
    "geom_block": "%geom MaxIter 50 Calc_Hess true end",
    "optimized_xyz": "S2_Retro/ts_guess.xyz",
    "output_file": "S2_Retro/preopt/ts_guess/...",
    "energy": -123.456,
    "imaginary_count": 1,
    "fell_back_to_raw": false
  },
  "intermediate": {
    "attempted": true,
    "accepted": true,
    "method_used": "Fast3C_Normal_PreOpt",
    "task_type": "opt",
    "route": "! r2SCAN-3c Opt tightSCF TightOpt CPCM(Acetone)",
    "optimized_xyz": "S2_Retro/intermediate.xyz",
    "energy": -124.000,
    "fell_back_to_raw": false
  }
}
```

## 11. checkpoint signature

`compute_step2_signature()` 中 `pre_opt` 应包含：

```python
"pre_opt": {
    "enabled": bool(...),
    "backend": "fast_3c_s3_style",
    "engine": pre_opt_cfg.get("engine"),
    "method": pre_opt_cfg.get("method"),
    "basis": pre_opt_cfg.get("basis"),
    "aux_basis": pre_opt_cfg.get("aux_basis"),
    "solvent": pre_opt_cfg.get("solvent"),
    "solvent_model": pre_opt_cfg.get("solvent_model"),
    "ts_guess": {
        "mode": ...,
        "route_task": ...,
        "max_cycles": ...,
        "calc_hess": ...,
        "freeze_forming_bonds": ...,
        "include_freq": ...,
        "require_one_imaginary": ...,
        "accept_zero_imaginary": ...,
    },
    "intermediate": {
        "mode": ...,
        "route_task": ...,
        "max_cycles": ...,
    },
}
```

这样切换 `r2SCAN-3c` / `B97-3c` 或 TS preopt 参数时，S2 cache 会失效。

## 12. 测试计划

### 12.1 单元测试

1. `test_s2_preopt_ts_uses_optts_not_opt`
   - mock ORCA 执行。
   - 断言 TS guess route 包含 `OptTS`。
   - 断言不等于普通 `Opt`。

2. `test_s2_preopt_keeps_fast3c_method`
   - 配置 `method: r2SCAN-3c`。
   - 断言 route 包含 `r2SCAN-3c`。
   - 配置 `method: B97-3c`。
   - 断言 route 包含 `B97-3c`。

3. `test_s2_preopt_intermediate_uses_normal_opt`
   - 断言 intermediate route 包含 `Opt`。
   - 断言 intermediate route 不包含 `OptTS`。

4. `test_s2_preopt_fallback_on_zero_imaginary`
   - TS preopt 返回 `imaginary_count=0`。
   - 默认 `accept_zero_imaginary=false`。
   - 断言回退 raw。

5. `test_s2_preopt_does_not_modify_s3_theory`
   - 运行 preopt 后检查原 config 中 `theory.optimization` 未被改写。

6. `test_step2_signature_changes_on_fast3c_method_change`
   - `r2SCAN-3c` 改为 `B97-3c`。
   - 断言 signature 改变。

### 12.2 输入文件快照测试

TS guess 应生成：

```text
! r2SCAN-3c OptTS Freq ...
```

或：

```text
! B97-3c OptTS Freq ...
```

不能出现：

```text
! r2SCAN-3c Opt ...
```

用于 TS guess。

### 12.3 端到端验证

样例：

```text
Output/la_zncl2_as_mgcl2/RXN_8ebf416f/branches/BR_MAJOR
```

验收：

1. `S2_Retro/preopt/ts_guess/*.inp` 包含 `OptTS`。
2. `S2_Retro/preopt/ts_guess/*.inp` 包含 `r2SCAN-3c` 或 `B97-3c`。
3. `S2_Retro/preopt/intermediate/*.inp` 包含 `Opt`，不含 `OptTS`。
4. `preopt_meta.json` 清楚记录 task_type、method、route、imaginary_count。
5. S3 继续消费 `S2_Retro/ts_guess.xyz`。
6. 若 TS preopt 坍塌到 0 虚频，S2 回退 `ts_guess_raw_s2.xyz`。

## 13. 分阶段实施

### P0：止血

1. 禁止 TS guess 使用当前 `constrained_optimize(task_type="opt")` 路径。
2. 若 `ts_guess.mode=optts` 但代码只能生成 `Opt`，直接 warning 并跳过 preopt，回退 raw。

### P1：ORCAInterface 通用任务入口

1. 新增 `optimize_with_task()`。
2. 支持 `opt`、`opt_freq`、`ts`、`ts_freq`。
3. 支持传入 `%geom MaxIter`。
4. 保持 composite_3c basis/aux 跳过逻辑。

### P2：QCTaskRunner fast-3c preopt API

1. 新增 `run_fast3c_ts_preopt_cycle()`。
2. 新增 `run_fast3c_normal_preopt_cycle()`。
3. 复用 S3 的 `_write_optimized_xyz()`、虚频统计和错误处理。

### P3：S2PreOptDriver 改造

1. 构造 fast-3c 临时 runner config。
2. 调用 QCTaskRunner preopt API。
3. 更新 `preopt_meta.json`。
4. 保持内嵌在 `run_step2()`。

### P4：测试和端到端验证

1. 单元测试通过。
2. import gate 通过。
3. 对 `RXN_8ebf416f/BR_MAJOR` 生成正确 fast-3c `OptTS` 输入。

## 14. 验收标准

必须满足：

1. 不新增 S2.5 数据流。
2. S2 对外 artifact 契约不变。
3. TS guess preopt 使用 fast-3c 方法。
4. TS guess preopt 使用 `OptTS` 或等价 TS 优化关键词。
5. TS guess preopt 不再使用普通 `Opt`。
6. intermediate preopt 可使用普通 `Opt`。
7. preopt 代码复用 S3-style runner/result/fallback 能力。
8. S3 最终 DFT 配置不被 preopt 修改。

最终推荐：

```text
S2 内嵌 preopt = fast-3c OptTS pre-pass
              + S3-style QCTaskRunner execution
              + raw fallback
              + no new pipeline stage
```
