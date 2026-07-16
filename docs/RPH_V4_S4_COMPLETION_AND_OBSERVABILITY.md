# RPH V4 S4 高精度计算与可观测性规范

**状态：** 已实施基础协议、S4 结构级日志与 WSL 终端查看器；branch-aware 输入扩展与高精度 IRC 为后续工作。  
**适用范围：** RPH V4 的 S4 `HighLevelEngine`，不适用于已经迁出的 `RPH_Postprocess` 特征提取与 ML 训练。

## 1. S4 的职责边界

S4 是所有候选结构的最终高精度几何与电子能层，不是“只重算一个 TS”的补丁步骤。它的设计目标是：

1. 为每个来自 S3 的结构建立独立的高精度 OPT、TS Freq（仅 TS）和 SP 记录；
2. 使高精度失败不抹掉 S3 已得到的低级别结构与能量；
3. 输出能被后处理/ML 安全消费的结构级状态、方法、溶剂与溯源；
4. 通过文件化事件流，使长时间 WSL 任务在终端断开、重新连接后仍可被准确观察。

S4 不承担以下职责：

- 不在 S4 重新识别反应机理、forming bonds 或立体化学；
- 不将某个分支自动判为“唯一正确路径”；
- 不以经验规则把真实 Lewis acid 的能量修正写入反应能垒；
- 不删除 S3 未收敛或 S4 未收敛的样本。

## 2. 固定计算协议

所有理论级别只由 `config/defaults.yaml` 给出，代码不得硬编码方法或溶剂。

| 结构角色 | 高精度 OPT | 高精度 Freq | 高精度 SP | 目的 |
|---|---|---|---|---|
| product / precursor / intermediate | ORCA M062X/def2-SVP, CPCM(acetone) | 默认不做 | ORCA wB97M-V/def2-TZVPP, CPCM(acetone) | 高精度极小值几何与电子能 |
| TS | ORCA M062X/def2-SVP `OptTS`, CPCM(acetone) | ORCA M062X/def2-SVP, CPCM(acetone), 独立 `Freq` 作业 | ORCA wB97M-V/def2-TZVPP, CPCM(acetone) | 高精度 TS 几何、虚频计数与电子能 |

S4 的 TS Freq 必须独立于 OptTS 运行，输入为 S4 优化后的几何。这避免将频率 Hessian 固定在优化前结构，也使失败能够单独记录。

### 溶剂约束

S4 的 ORCA 优化、频率与单点作业均要求显式 CPCM。实际输入路由必须包含：

```text
CPCM(acetone)
```

仅在 YAML 中写 `solvent: acetone` 而没有将其渲染进 ORCA route，等价于气相计算，不能接受。S4 不默认使用 SMD；若未来要变更模型，必须同时修改 YAML、计算签名、基准测试与 manifest 方法记录。

## 3. S3 → S4 结构选择与溯源

S4 对每个结构使用以下几何优先级：

```text
S3 opt_xyz  →  S3 input_xyz fallback  →  fail-fast（无 XYZ）
```

`S3_fallback` 是显式降级状态，不代表高精度优化成功。每个 S4 manifest 条目必须保留：

```json
{
  "id": "ts_peb_peak",
  "kind": "ts",
  "input_source": "S3_OPT",
  "source_s3": {
    "structure_status": "complete",
    "opt_status": "complete",
    "sp_status": "complete",
    "frequency_status": "complete",
    "ts_frequency_valid": true,
    "manifest": "/tmp/rph4_rx_1/S3_LowLevel/manifest.json"
  }
}
```

后续 branch 重构完成后，条目还必须包含 `reaction_id`、`condition_id`、`branch_id`、`pathway_id`、`role`、`parent_structure_id` 与原子映射空间。S4 不得通过目录名猜测这些信息。

## 4. S4 状态语义与 ML 使用规则

| `status` | 定义 | `usable_for_ml` |
|---|---|---|
| `complete` | OPT 和 SP 成功；若为 TS，Freq 也通过单一显著虚频门槛 | `true` |
| `ts_frequency_unverified` | OPT、SP 成功，但 TS Freq 失败或虚频数不正确 | `true`，但必须带 TS 质量标签 |
| `opt_failed_sp_complete` | OPT 失败，但从输入/fallback 几何完成 SP | `true`，仅作为降级电子结构样本 |
| `degraded` | OPT 或 SP 中至少一个失败，且未满足上一行 | `false` |
| `failed` | 作业级异常、缺输入或未产生可用 SP | `false` |

`usable_for_ml=true` 的含义只是“存在可用的该级别 SP”，并不表示它是已验证反应路径。下游训练必须同时使用：

- `s4_status`；
- `input_source`；
- `source_s3.ts_frequency_valid`；
- S4 `ts_frequency_valid`；
- 后续 branch/condition 身份。

当前 S4 Freq 的自动门槛是“恰好一个 `<= -50 cm⁻¹` 的虚频”。它验证的是虚频计数，**尚不等价于反应坐标位移验证或 IRC 验证**。这两个字段必须分别建模，不能用一个 `ts_valid` 布尔值混在一起。

## 5. 能量比较规则

只有以下身份相同的结构能量允许直接比较：

```text
同一 reaction_id
+ 同一 condition_id
+ 同一 branch_id / pathway_id
+ 同一 charge / multiplicity
+ 同一方法、基组、溶剂模型
+ 相同能量定义（SP electronic energy 或含热校正自由能）
```

S4 的 wB97M-V SP 是电子能。除非显式引入同级别热校正，不得把它命名为 `ΔG‡`。S3 的热化学、S4 的电子能和实验选择性标签应作为不同字段进入后处理。

## 6. 输出契约

每次真正执行 S4 后，`S4_HighLevel/` 必须包含：

```text
S4_HighLevel/
├── status.json       # 当前完整状态快照；原子替换写入
├── events.jsonl      # 追加式事件流；一行一个 JSON 事件
├── s4.log            # 人类可读的同源事件日志
├── manifest.json     # S4 最终结构结果与方法快照
└── <structure_id>/
    ├── opt/
    ├── freq/         # 仅 TS，或在配置明确要求时
    └── sp/
```

运行根目录额外产生 `rph_v4.log`；可用 CLI 的 `--log-file` 覆盖。

### `status.json`

`status.json` 不是最终结果文件，而是 WSL 查看器读取的实时快照。必须有：

- 总体 `status`：`running`、`completed` 或 `completed_with_failures`；
- `started_at`、`updated_at`、`finished_at`，统一 ISO-8601 UTC；
- `summary.total/pending/running/finished/failed/usable_for_ml`；
- 每个结构的 `current_task`、各任务状态、错误、SP 能量与 TS Freq 结果。

写入采用临时文件后原子替换，查看器不会读到半截 JSON。

### `events.jsonl`

事件至少覆盖：

```text
stage_started / stage_finished
structure_started / structure_finished / structure_failed
optimization_started / optimization_finished
frequency_started / frequency_finished
frequency_skipped
single_point_started / single_point_finished
```

每个任务事件记录 `structure_id`、engine、method、solvent、solvent_model、完成状态、输出文件、能量和错误。该文件是排查 ORCA 失败、重建时间线和后续 UI 的唯一增量数据源。

## 7. WSL 实时可视化

第一终端启动计算：

```bash
python bin/rph_run \
  --csv "$RPH_ROOT/data/reaxys_cleaned.csv" \
  --rx-id 1 \
  --output /tmp/rph4_rx_1 \
  --config "$RPH_ROOT/config/defaults.yaml" \
  --stop-after s4
```

第二终端查看动态 S4 面板：

```bash
python bin/rph_watch \
  --output /tmp/rph4_rx_1 \
  --watch \
  --interval 3
```

该面板显示总体完成数、失败数、当前结构、当前 QC 子任务和 OPT/Freq/SP 状态。它不读取正在写入的 ORCA 原始输出，因此不会与计算进程争抢文件。

需要原始事件时间线时：

```bash
tail -F /tmp/rph4_rx_1/S4_HighLevel/s4.log
tail -F /tmp/rph4_rx_1/rph_v4.log
```

任务结束后，可打印一次静态摘要：

```bash
python bin/rph_watch --output /tmp/rph4_rx_1 --no-color
```

## 8. 与后续 branch 架构的衔接

当前实现的 S4 日志协议已经按“多个独立结构”设计，但现有 V4 编排器尚未恢复完整 DR/condition/precursor branch 拓扑。完成 branch 重构时，S4 输入应扩展为：

```text
RXN/<reaction_id>/
  branches/<branch_id>/S4_HighLevel/
  conditions/<condition_id>/branches/<branch_id>/S4_HighLevel/
```

每个 branch 独立产生自己的 `status.json`、`events.jsonl` 和 checkpoint。顶层反应状态只汇总子 branch，不能将不同 branch 的能量、失败或缓存混写到同一个结构条目。

每个 branch 的 S4 最低任务数为：

```text
2 × (N_product + N_precursor + N_intermediate)
+ 3 × N_TS
```

其中 TS 的三项为 OptTS、Freq、SP。这个公式将取代当前扁平 V4 中遗漏 precursor 的任务清单。

## 9. 验收标准

一次 S4 代码变更只有同时满足以下条件才算完成：

1. ORCA S4 输入中可见 `CPCM(acetone)`；
2. ORCA S4 SP manifest 记录 CPCM；
3. TS 在高精度 OptTS 后确实出现独立 `freq/` 作业和频率数组；
4. 单个结构失败不会阻止下一个结构继续执行；
5. `status.json` 可在作业运行期间持续读取；
6. `events.jsonl` 包含任务开始与结束事件；
7. `bin/rph_watch` 可在无 Rich、无 jq 的原生 WSL 环境工作；
8. 配置或 S3 输入发生变化时，S4 checkpoint 签名失效并重跑；
9. S4 manifest 保留 S3 溯源与降级输入来源；
10. 回归测试覆盖 ORCA CPCM route、ORCA Freq、S4 日志和查看器。

## 10. 后续增强（不在当前自动化中隐式开启）

- 高精度 TS 正常模位移与 forming-bond 投影验证；
- 仅对通过 S3/S4 门槛的 TS 启动高精度 IRC；
- 分支级并发调度与资源预算；
- 将 `events.jsonl` 接入真正的 GUI/Web UI；
- 由 RPH_Postprocess 消费 S4 manifest，生成 branch-aware 特征与 ML 数据集。

这些功能必须建立在本文件的 manifest、状态与溯源契约之上，不能回到解析临时 stdout 或猜测目录结构的旧实现。
