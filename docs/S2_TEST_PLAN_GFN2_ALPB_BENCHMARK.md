# S2 GFN2-xTB/ALPB Benchmark Test Plan

**目标目录**: `RPH_Test_Results/s2_gfn2_benchmark_alpb_v2`（20 个反应，rx_id_1..20）
**流水线范围**: `--stop-after s2`（S2 PEB 扫描，不进入 S3/S4）
**S2 方法**: ORCA GFN2-xTB relaxed scan（ALPB 溶剂）→ B97-3c SP 全帧精修 → 种子选择
**状态日期**: 2026-08-04

> **当前模式（2026-08-04 更新）**：GFN2-only 鲁棒性测试优先——`step2.rescue.enabled: false`，
> 拓扑漂移不再触发 B97-3c rescue，而是容忍（保留 knee 移位的 GFN2 主选择，记录
> `orca_gfn2_topology_drift_tolerated_no_rescue` 降级标记）。全量 20 反应先验证 GFN2 独立鲁棒性，
> 通过后再评估是否恢复 rescue（Phase B 之后）。

---

## 1. 背景与目标

三个历史故障已修复并实测闭环：

| # | 故障 | 修复 | 验证 |
|---|------|------|------|
| 1 | GFN2-xTB 输入带 `CPCM(Acetone)`，ORCA abort (rc=25) | `_render_route` 支持 ALPB；`orca_gfn2_scan` 配置 `solvent_model: ALPB`；semiempirical_xtb 自动翻译 CPCM→ALPB | rx_id_1 扫描输入为 `ALPB(Acetone)` ✓ |
| 2 | ORCA 5.0.4 找不到 xtb（只查自身目录，不查 PATH） | `ln -s /opt/software/xtb/bin/xtb /opt/software/orca/xtb`；`run_surface_scan` 预检 + `XTB_BINARY_UNAVAILABLE` 分类 | 17 帧全部产出 ✓ |
| 3 | 孤儿 ORCA 进程抢 16 核 | 清理（pkill 遗留 worker） | 扫描 6m19s（争抢时 20min+）✓ |

**性能基线**（ORCA 5.0.4，rx_id_1 product_major，41 原子）：
- 17 扫描步 × ~14 几何循环 = **238 次串行 xtb 调用**
- xtb 单次 SCF 实测 **0.028-0.05s**（1/4/8/16 线程零加速 → 16 核对该分子无用）
- ORCA 侧每调用开销 **~1.55s**（fork 子进程 + xcontrol 写读握手）→ 占墙钟 **94%**
- 总墙钟 6m32s：真实计算仅 22.3s（xtb 10.5s + ORCA 几何 11.8s）

**本方案的验收对象**：S2 全流程正确性（20 反应）+ 性能假设验证（ORCA 6.1.0 是否大幅降低接口开销）。

---

## 2. 前置门禁（Phase 0，全部必须通过）

| 门禁 | 检查 | 命令 |
|------|------|------|
| G0.1 | 代码编译 + 导入门禁 | `python -m py_compile rph_core/utils/orca_interface.py rph_core/utils/orca_failure_classifier.py rph_core/steps/step2_retro/relaxed_scan_rescue.py`；`python scripts/ci/check_imports.py rph_core` |
| G0.2 | 单元测试 | `pytest -q tests/test_s2_orca_gfn2_workflow.py tests/test_v4_orca_convergence.py tests/test_v4_protocol_contract.py tests/test_v4_checkpoint.py tests/test_v4_stage_calculator.py` |
| G0.3 | xtb 位于 ORCA 目录 | `ls /opt/software/orca/xtb`（符号链接存在） |
| G0.4 | 无孤儿 QC 进程 | `ps -eo pid,cmd \| grep -E "orca\|xtb" \| grep -v grep`（除当前作业外为空） |
| G0.5 | config 冻结 | `config/defaults.yaml` 中 `step2.orca_gfn2_scan.relaxed_scan` 含 `solvent: acetone` + `solvent_model: ALPB` |

> G0.5 注意：首次 `bin/rph_run` 会冻结 `run.config.json`。修改 config 后需确认目标目录的 `run.config.json` 同步（必要时删除重建，或依赖 `--refresh-stages` 的 hash 校验机制报差异）。

---

## 3. 正确性测试（Phase A：全量 20 反应，GFN2-only 模式）

### 3.1 执行

```bash
CSV=data/reaxys_cleaned.csv
NEW=RPH_Test_Results/s2_gfn2_benchmark_alpb_v2
# 4 核/反应 × N 并行（8 = 32 核满载）；rescue 已在 config 禁用
bash scripts/run_s2_parallel_4core.sh "$CSV" "$NEW" 8
# 脚本自动产出 $NEW/SUMMARY.csv（每变体：s2_state / 种子索引 / 排除帧 / 拓扑决策）
```

### 3.2 每反应验收门禁（C1-C7）

对每个 `S2_PEB/<variant>/manifest.json`：

| 门禁 | 内容 | 通过标准 |
|------|------|---------|
| C1 | 扫描 route | 输入的 `.inp` 含 `ALPB(` 且无 `CPCM(`；`%pal nprocs 4` |
| C2 | GFN2 扫描 | `scan_method=GFN2-xTB`，帧数 ≥ 3，`s2_state ∈ {gfn2_seeded, unresolved}`（GFN2-only 模式下 `rescue_seeded` 不应出现） |
| C3 | B97-3c SP 全覆盖 | `energy_refinement_method=B97-3c`，精修帧数 = 扫描帧数 |
| C4 | 种子产物 | `ts_guess.xyz` + `intermediate.xyz` 均存在且非空（≥ 3 原子） |
| C5 | manifest | `status=COMPLETE`，`mapping_status=verified`，变体 manifest 齐全 |
| C6 | 拓扑漂移处理 | 若有 `Topology drift` 警告：`excluded_frames` 明确记录；**不得创建 rescue 目录**；`topology_state` 与 `rescue_decision` 可解释（tail_distorted_after_knee / rescue_*）；种子索引与漂移帧不相交 |
| C7 | 可诊断性 | 若失败（unresolved），日志含明确分类（`XTB_BINARY_UNAVAILABLE`/`orca_gfn2_scan_or_sp_incomplete` 等），非笼统报错 |

### 3.3 GFN2 鲁棒性验收标准（本阶段核心）

基于 `SUMMARY.csv` 全量统计：

- **硬门槛**：`gfn2_seeded` 率 ≥ 85%（≥ 17/20 反应至少一个变体成功选种）
- **漂移容忍率**：记录 `excluded_frames` 非空但仍 `gfn2_seeded` 的变体数（应为多数漂移案例）
- **unresolved 上限**：≤ 15%（>3 反应完全失败即失败批次），且每个必须可诊断（C7）
- **rescue 禁用确认**：`rescue_seeded` 计数必须为 0；无任何反应创建 rescue 目录
- 每个 unresolved 单独诊断：GFN2 扫描失败（真失败）vs 拓扑全排除导致无种（GFN2 能力边界）
- 输出：`SUMMARY.csv` + 漂移分布统计（漂移帧索引直方图，验证是否集中在尾段 14-16）

---

## 4. 性能测试（Phase B：验证 16 核与 ORCA 6.1.0 假设）

### B1. nproc A/B（证明 16 核不是瓶颈，释放核数）

- 在 `step2.orca_gfn2_scan.relaxed_scan` 增加 `nproc: 4`（需代码支持：`relaxed_scan_rescue.py` 将 `scan_cfg["nproc"]` 传入 spec）
- 同一反应（rx_id_2）分别以 nproc=16 / 4 跑 GFN2 扫描
- **通过标准**：墙钟差 ≤ 15%（预期基本持平），且 nproc=4 时系统其余核可并行跑别的反应

### B2. ORCA 6.1.0 A/B（核心假设验证）

前置：
1. 安装 ORCA 6.1.0 到独立目录（如 `/opt/software/orca6`），**不动 5.0.4**
2. 确认 6.1.0 自带 `otool_xtb`（ORCA 6.0.1+ 捆绑）或按官方要求放置 xtb（6.1 手册要求 xtb ≥ 6.7.1）
3. `config` 切换 `executables.orca.path` → ORCA 6.1.0，`ld_library_path` 同步

同一反应（rx_id_2）跑 5.0.4 vs 6.1.0 对比：

| 指标 | 采集方式 |
|------|---------|
| 总扫描墙钟 | `TOTAL RUN TIME`（.out） |
| 每调用开销 | `(TOTAL - Sum of individual times) / SCF 调用数` |
| 能量一致性 | 两版本 `scan_profile.json` 的 `energies_hartree`（B97-3c 精修后）逐帧 RMSD |
| 种子一致性 | `ts_seed_index` / `int_seed_index` 是否相同 |

- **通过标准（性能）**：6.1.0 墙钟 ≤ 5.0.4 的 1/3（若接口确实进程内化，预期 5-10x）
- **通过标准（科学一致性）**：B97-3c 精修能量 RMSD < 0.05 kcal/mol（同一 xtb 内核下应几乎为 0）；种子索引一致；`s2_state` 一致
- **若 B2 不通过**（6.1.0 开销未降）：性能改进转入"原生 xtb 扫描"方案（`xtb_runner.run_scan()`，进程内约束优化）——此方案改变科学路径，需先做能量剖面等价性验证再采纳

---

## 5. 鲁棒性测试（Phase C：故障注入与降级路径）

| 测试 | 方法 | 预期 |
|------|------|------|
| F1 | 拓扑漂移 → GFN2-only 容忍 | 利用漂移反应（如 rx_id_1/2/3 的 frames 14-16）：确认不创建 rescue 目录、`s2_state=gfn2_seeded`、`degraded_reasons` 含 `orca_gfn2_topology_drift_tolerated_no_rescue`、`topology_state` 可解释 |
| F2 | xtb 缺失预检 | 临时移走 `/opt/software/orca/xtb`，跑单反应，确认 `run_surface_scan` 预检立即失败且错误含 `ln -s` 提示（不白跑 ORCA） |
| F3 | 断点恢复 | 中途 Ctrl-C 后重跑同一命令，确认 `pipeline.state` 恢复、已完变体 reuse |
| F4 | 资源隔离 | 两个反应并行跑，确认调度不超卖、无孤儿进程残留（对照 G0.4） |

---

## 6. S2→S3 交接验证（Phase D）

S2 完成后随机抽 2 个反应（如 rx_id_1、rx_id_3）执行：

```bash
bin/rph_run --config config/defaults.yaml --csv "$CSV" --rx-id 1 \
  --output "$NEW/rx_id_1" --stop-after s3
```

- **H1**：`ts_guess.xyz` / `intermediate.xyz` 被 S3 正确消费（结构可读、元素数匹配）
- **H2**：S3 中 `forming_bonds`（XYZ 索引）与 S0 manifest 一致
- **H3**：S3 能启动 `OptTS`（TS 种子）与 `MIN`（INT 种子）作业——不要求收敛，只要求正确启动
- 该阶段通过 = S2 产物可交接性证明，同时暴露种子质量问题（如 TS 种子需大幅优化）

---

## 7. 汇总报告输出

每个反应的最终记录写入 `RPH_Test_Results/s2_gfn2_benchmark_alpb_v2/SUMMARY.csv`：

```csv
rx_id, variant, s2_state, selection_source, selection_method, ts_seed_index, int_seed_index,
scan_elapsed_s, sp_frames, sp_covered, excluded_frames, confidence, degraded, notes
```

汇总后按 3.3 的标准判定批次通过/失败，并将对比数据（5.0.4 vs 6.1.0）附在 `docs/` 下的性能对比小节。

---

## 8. 执行顺序与工时预估（单核 16 线程机器）

| 阶段 | 内容 | 预计耗时 |
|------|------|---------|
| 0 | 门禁 G0.1-G0.5 | 15 min |
| A | 19 反应批量（每反应 ~6-10 min，串行） | 2-3 h |
| B1 | nproc A/B | 30 min |
| B2 | ORCA 6.1.0 安装 + A/B | 1-2 h |
| C | 故障注入 4 项 | 1 h |
| D | S3 交接 2 反应 | 1-2 h（S3 实际计算） |
| 汇总 | 报告 | 30 min |

**建议并行**：Phase A 中每 4 个反应一组（16 核 ÷ 4 核/作业，若 B1 验证 nproc=4 无损则全组并行），总批次时间可压到 ~1 h。

---

## 9. 风险与已知限制

1. **ORCA 6.1.0 性能提升未验证**——B2 是本方案的关键实验，若 6.1.0 仍为每调用子进程，开销改善有限（librarian 结论见方案配套分析）
2. **ALPB vs CPCM 的科学差异**——GFN2 扫描用 ALPB 而 S3/S4 用 CPCM，溶剂模型不一致，需确认这是可接受的基准设定（ALPB 是 GFN-xTB 唯一可用选项，属 ORCA 硬约束，非可选项）
3. **B97-3c SP 精修开销**：17 帧 × 每帧 1-3 min（16 并行）≈ 2-4 min/反应——这是与 GFN2 扫描同量级的第二耗时项，后续若需加速可评估 SP 缓存/降帧
4. **run.config.json 冻结机制**：切换 ORCA 6.1.0 后 config 变化会触发 checkpoint 差异，需按 README 的 `--resume-policy`/`--recompute-from` 语义处理
