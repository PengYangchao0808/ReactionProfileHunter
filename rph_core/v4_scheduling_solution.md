# S1 并行调度方案：吞吐优先与 Profiled Final Wave

## 1. 已确认的问题

`rx_id_1_test3` 的 `product_minor_001` 提供了直接证据：同一批 B97-3c
单点中，1 核作业约需 1 分 50 秒至 2 分 06 秒，而 2 核作业约需 3 分 08
秒至 4 分 56 秒。SCF 循环数均为 13–15，因此主要问题不是收敛，而是 ORCA
多 rank 负扩展、作业并发和资源竞争。

旧调度器在启动前静态执行：

```text
bulk = records[:-8]
drain = records[-8:]
```

对于 12 个构象，这会在批次开始时形成 `4×1 + 6×2`，并让最后两个任务
排队。它不是真正的运行期尾部调度。

## 2. 当前生产默认值

S1 variant 仍顺序执行，每个 variant 可以使用完整资源预算：

```text
CREST/GFN2              configured cores
B97-3c bulk             up to 16 jobs × 1 rank
xTB mRRHO bulk          up to 16 jobs × 1 thread
```

允许核心空闲。候选数为 12 时，默认运行 `12×1`，不会为了填满 16 核而将
作业扩为 2 核。

`queue_drain` 和 `adaptive_tail` 已废弃。即使旧配置仍包含并启用这些字段，
运行时也只给出警告并使用安全的 throughput-first 模式。

## 3. Profiled final wave

调度分为两个严格隔离的阶段：

```text
throughput bulk
    所有任务固定 cores_per_job
    持续填充 parallel_jobs
    |
    | barrier：bulk 全部结束
    v
profiled final wave（可选）
    只接受固定 ensemble benchmark 得到的精确布局
```

final wave 不能与 bulk 混跑。它按 `candidate_count % parallel_jobs` 确定最后
一波任务数；若余数为零，则最后一波为一个完整并行波次。

只有同时满足以下条件才启用：

- `scheduling.final_wave.enabled=true`；
- `profile_layouts` 中存在完全匹配的任务数；
- 布局长度等于任务数；
- 总核数不超过全局预算；
- 单任务核数不超过配置上限；
- 布局满足并发作业数限制。

示例（仅在 benchmark 证明后配置）：

```yaml
b97_3c:
  parallel_jobs: 16
  cores_per_job: 1
  scheduling:
    policy: throughput_first
    final_wave:
      enabled: true
      max_cores_per_job: 4
      profile_layouts:
        "6": [3, 3, 3, 3, 2, 2]
```

没有精确 profile、profile 无效或 profile 与单核布局相同，都会回退到连续
单核吞吐模式，不设置 barrier。

## 4. ORCA 内存预算

旧逻辑只按单个 ORCA 作业的 rank 数计算 `%maxcore`，没有考虑同时运行的作业
数量。现在每个调度波次统一计算：

```text
usable_memory = total_memory × orca_maxcore_safety
rank_memory = usable_memory - job_count × job_overhead
maxcore_per_rank = rank_memory / total_active_ranks
```

若每 rank 的 `maxcore` 低于 `min_maxcore_mb_per_rank`，调度器先降低 ORCA
并发作业数；profiled final wave 若仍不能满足预算则拒绝运行。

默认配置：

```yaml
scheduling:
  memory:
    job_overhead_mb: 256
    min_maxcore_mb_per_rank: 1000
    maxcore_cap_mb_per_rank: 4000
```

`%maxcore` 是上限而不是实际 RSS；真实峰值内存仍需由 WSL benchmark 记录。

## 5. 隐式线程隔离

ORCA 的 rank 数由 `%pal nprocs` 控制。每个 ORCA 进程启动时显式设置：

```text
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
OMP_MAX_ACTIVE_LEVELS=1
OMP_THREAD_LIMIT=1
```

这样底层 BLAS/OpenMP 不会在调度器统计之外创建额外线程。每 rank 线程数可通过
`resources.orca_math_threads_per_rank` 配置，但默认值为 1。

## 6. 可观测性

每个 B97-3c 候选记录：

- scheduler layout；
- submit/start/finish 时间；
- queue wait 和 wall time；
- `nprocs`；
- `%maxcore`；
- ORCA `TOTAL RUN TIME`；
- SCF 循环数；
- 输出文件路径。

mRRHO 记录同样的调度、队列和 wall-time 字段。UI batch 事件包含调度策略与
final-wave 布局。

## 7. Benchmark 门禁

final-wave profile 必须来自固定 XYZ ensemble，不能比较重新运行 CREST 后得到的
不同构象集合。对 4、6、8、12、16、32 和大 ensemble 分别比较 `N×1`、`8×2`、
`4×4` 及候选混合布局，至少记录：

- batch wall time；
- total core-hours；
- P50/P95 单任务时间；
- peak RSS 和 swap；
- CPU utilization 和 I/O wait；
- SCF/SCC failure rate。

在 WSL T0–T6 门禁通过前不运行大型 benchmark。未通过 benchmark 的布局不得写入
生产 profile。

## 8. 科学语义与 checkpoint

调度不会改变候选集合、B97-3c 能量、mRRHO 修正或 ensemble partition function。
`parallel_jobs`、核数、`scheduling`、内存调度和 ORCA rank 内线程数均不进入 S1
科学签名，因此仅调整性能布局不会无意义地使科学 checkpoint 失效。
